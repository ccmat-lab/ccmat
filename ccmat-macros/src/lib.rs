/*
*
* The **math related** macros defined in the crate are major for internal usage.
* The macros that for users will be re-exported by `ccmat`.
*/

use proc_macro::TokenStream;
use quote::quote;
use syn::parse::{Parse, ParseStream};
use syn::spanned::Spanned;
use syn::{BinOp, ExprBinary, ExprLit, Lit};
use syn::{Expr, Result, Token};

#[rustfmt::skip]
struct MatrixInput {
    mat: [[Expr; 3]; 3],
}

fn detect_integer_division(expr: &Expr) -> Result<bool> {
    match expr {
        Expr::Binary(ExprBinary {
            left, op, right, ..
        }) => {
            // If op is `/`, check if both sides look integer-ish
            if matches!(op, BinOp::Div(_)) && (is_integer_expr(left)? || is_integer_expr(right)?) {
                return Ok(true);
            }
            // Recurse down left/right anyway
            Ok(detect_integer_division(left)? || detect_integer_division(right)?)
        }

        // Parentheses
        Expr::Paren(p) => detect_integer_division(&p.expr),

        // Unary ops: `-1/2`
        Expr::Unary(u) => detect_integer_division(&u.expr),

        // Function call arguments, include Path Call e.g f64::cos(1/2)
        Expr::Call(call) => {
            let mut iflag = false;
            for arg in &call.args {
                if detect_integer_division(arg)? {
                    iflag = true;
                }
            }
            Ok(iflag)
        }

        // // Method call args: foo.bar(1/2)
        // Expr::MethodCall(call) => {
        //     detect_integer_division(&call.receiver) || call.args.iter().any(detect_integer_division)
        // }

        // Other expressions — no division here
        _ => Ok(false),
    }
}

fn is_integer_expr(expr: &Expr) -> Result<bool> {
    match expr {
        Expr::Path(_) => Err(syn::Error::new(
            expr.span(),
            "unable to know if it is a float when expanding macro",
        )),
        Expr::Lit(ExprLit {
            lit: Lit::Int(_), ..
        }) => Ok(true),
        Expr::Paren(p) => is_integer_expr(&p.expr),
        Expr::Unary(u) => is_integer_expr(&u.expr), // handle -1
        _ => Ok(false),
    }
}

impl Parse for MatrixInput {
    fn parse(input: ParseStream) -> Result<Self> {
        let mut rows = Vec::new();
        let mut current_row = Vec::new();
        let span = input.span();

        while !input.is_empty() {
            if let Ok(expr) = input.parse::<Expr>() {
                // detect integer devide as an error e.g raise on 1/2
                if detect_integer_division(&expr)? {
                    return Err(syn::Error::new(
                        expr.span(),
                        "integer division detected in matrix literal; e.g use `1./2.` instead of `1/2`",
                    ));
                }

                if current_row.len() > 2 {
                    return Err(syn::Error::new(expr.span(), "expect 3 items per row"));
                }
                current_row.push(expr);
            }

            if input.peek(Token!(;)) {
                input.parse::<Token!(;)>()?;
                rows.push(std::mem::take(&mut current_row));
            } else if input.peek(Token!(,)) {
                // allow ',' commas optionally
                input.parse::<Token!(,)>()?;
            } else {
                continue;
            }
        }

        // after the last
        if !current_row.is_empty() {
            rows.push(current_row)
        }

        if rows.len() != 3 {
            return Err(syn::Error::new(span, "expect 3 rows for a 3x3 matrix"));
        }

        let mat = [
            [rows[0][0].clone(), rows[0][1].clone(), rows[0][2].clone()],
            [rows[1][0].clone(), rows[1][1].clone(), rows[1][2].clone()],
            [rows[2][0].clone(), rows[2][1].clone(), rows[2][2].clone()],
        ];
        Ok(MatrixInput { mat })
    }
}

#[proc_macro]
pub fn matrix_3x3(tokens: TokenStream) -> TokenStream {
    let tokens: proc_macro2::TokenStream = tokens.into();
    let mat_input = syn::parse2::<MatrixInput>(tokens);

    let mat_input = match mat_input {
        Ok(input) => input,
        Err(err) => return err.to_compile_error().into(),
    };

    let MatrixInput { mat } = mat_input;

    let row0 = mat[0].iter().map(|x| quote!(f64::from(#x)));
    let row1 = mat[1].iter().map(|x| quote!(f64::from(#x)));
    let row2 = mat[2].iter().map(|x| quote!(f64::from(#x)));

    let expand = quote! {{
        let mat: [[f64; 3]; 3] = [
            [#(#row0,)*],
            [#(#row1,)*],
            [#(#row2,)*],
        ];
        mat
    }};

    expand.into()
}
