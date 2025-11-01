use ccmat_macros::matrix_3x3;

fn main() {
    let a = 2;
    let b = 4;
    let mat = matrix_3x3![
        a/b 2 3;
        4 5 6.1;
        7 8 9;
    ];
}
