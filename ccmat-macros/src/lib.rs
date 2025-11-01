/*
*
* The **math related** macros defined in the crate are major for internal usage.
* The macros that for users will be re-exported by `ccmat`.
*/

// struct MatrixInput {
//
// }
//
// #[proc_macro]
// pub fn matrix_3x3(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
//     let input = proc_macro2::TokenStream::from(input);
//
//     let output: proc_macro2::TokenStream = {
//         dbg!(input)
//     };
//     proc_macro::TokenStream::from(output)
// }

use proc_macro::TokenStream;
use syn::parse::{Parse, ParseStream};
use syn::{parse_macro_input, parse_quote, Expr, Result, Token};

#[rustfmt::skip]
struct MatrixInput {
    mat: [[Expr; 3]; 3],
}

impl Parse for MatrixInput {
    fn parse(input: ParseStream) -> Result<Self> {
        let mut rows = Vec::new();
        let mut current_row = Vec::new();

        while !input.is_empty() {
            if let Ok(expr) = input.parse::<Expr>() {
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
            row.push()
        }

        let mat = [
            [rows[0][0], rows[0][1], rows[0][2]],
            [rows[0][0], rows[0][1], rows[0][2]],
            [rows[0][0], rows[0][1], rows[0][2]],
        ];
        Ok(mat)
    }
}

#[proc_macro]
pub fn matrix_3x3(tokens: TokenStream) -> TokenStream {
    let input = parse_macro_input!(tokens as MatrixInput);


    todo!()
}
