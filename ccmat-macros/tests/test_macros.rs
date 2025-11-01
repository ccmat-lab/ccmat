use ccmat_macros::matrix_3x3;

#[test]
fn matrix_3x3() {
    let mat = matrix_3x3![
        1 2 3;
        4 5 6.1;
        7 8 9;
    ];

    assert_eq!(mat[0][0], 1.0);
    assert_eq!(mat[1][2], 6.1);

    let mat = matrix_3x3![
        1 + 2 + 1 2 3;  // althrough ugly and not recommemded
        4 5 6.1;
        7 8 9;
    ];

    assert_eq!(mat[0][0], 4.0);
    assert_eq!(mat[1][2], 6.1);

    let mat = matrix_3x3![
        1 + 2 + 1, 2, 3; // better
        4, 5, 6.1;
        7, 8, 9;
    ];

    assert_eq!(mat[0][0], 4.0);
    assert_eq!(mat[1][2], 6.1);
}

#[test]
fn matrix_3x3_frac_pass() {
    let a = 2;
    let b = 4;
    let mat = matrix_3x3![
        a/b 2 3;
        4 5 6.1;
        7 8 9;
    ];
    assert_eq!(mat[0][0], 0.5);
    assert_eq!(mat[1][2], 6.1);
}

// #[test]
// fn matrix_3x3_frac_error() {
//     let mat = matrix_3x3![
//         2/4 2 3;
//         4 5 6.1;
//         7 8 9;
//     ];
// }

// #[test]
// fn matrix_3x3_f64_func() {
//     let tmatrix = matrix_3x3![
//         sqrt(3)/2  -1/2     0;
//         1/2   sqrt(3)/2     0;
//         0             0     1;
//     ];
//
//     assert_eq!(mat[0][0], 0.5);
//     assert_eq!(mat[1][2], 6.1);
// }

