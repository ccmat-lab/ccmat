/*
*
* The **math related** macros defined in the crate are major for internal usage.
* The macros that for users will be re-exported by `ccmat`.
*/

use proc_macro::TokenStream;
use quote::quote;
use syn::parse::{Parse, ParseStream};
use syn::spanned::Spanned;
use syn::{Error, Expr, Result, Token};

#[rustfmt::skip]
struct MatrixInput {
    mat: [[Expr; 3]; 3],
}

impl Parse for MatrixInput {
    fn parse(input: ParseStream) -> Result<Self> {
        let mut rows = Vec::new();
        let mut current_row = Vec::new();
        let span = input.span();

        while !input.is_empty() {
            if let Ok(expr) = input.parse::<Expr>() {
                if current_row.len() > 2 {
                    return Err(Error::new(expr.span(), "expect 3 items per row"));
                }
                current_row.push(expr)
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
            return Err(Error::new(span, "expect 3 rows for a 3x3 matrix"));
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
    let MatrixInput { mat } = syn::parse2(tokens).expect("failed to parse input matrix");
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
