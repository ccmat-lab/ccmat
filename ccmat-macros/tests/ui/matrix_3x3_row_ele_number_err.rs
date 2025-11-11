use ccmat_macros::matrix_3x3;

fn main() {
    let tmatrix = matrix_3x3![
        sqrt(3.)/2.    -1./2.     0;
        1./2.     sqrt(3.)/2.     0;
        0                  0      1;
    ];
}
