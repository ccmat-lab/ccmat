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
        1, 2, 3;
        4, 5, 6.1;
        7, 8, 9;
    ];

    assert_eq!(mat[0][0], 1.0);
    assert_eq!(mat[1][2], 6.1);

    // let tmatrix = matrix_3x3![
    //     sqrt(3)/2  -1/2     0;
    //     1/2   sqrt(3)/2     0;
    //     0             0     1;
    // ];
}
