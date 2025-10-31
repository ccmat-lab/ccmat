use std::ops::{Add, Deref, Index, Mul};

use crate::structure::{Angstrom, InvAngstrom};

#[derive(Default)]
pub struct Matrix3(pub [[f64; 3]; 3]);

impl Index<usize> for Matrix3 {
    type Output = [f64; 3];

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

pub type TransformationMatrix = Matrix3;
pub type RotationMatrix = Matrix3;

#[macro_export]
macro_rules! matrix_3x3 {
    ( $( $( $x:expr )+ );+ $(;)? ) => {{
        let inner = [
            $(
                [ $( f64::from($x)),+ ],
            )+
        ];
        $crate::math::Matrix3(inner)
    }};
}

// TODO: add handy idx accessing method
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vector3<T>(pub [T; 3]);

impl<T> Add<Vector3<T>> for Vector3<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, rhs: Vector3<T>) -> Self::Output {
        Vector3([
            (*self)[0] + (*rhs)[0],
            (*self)[0] + (*rhs)[0],
            (*self)[0] + (*rhs)[0],
        ])
    }
}

impl<T> Deref for Vector3<T> {
    type Target = [T; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> Index<usize> for Vector3<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Mul<Vector3<f64>> for f64 {
    type Output = Vector3<f64>;

    /// Scalar multiply for a `Vector3<f64>`
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmat_core::math::Vector3;
    ///
    /// let vec3 = Vector3::<f64>([2.0, 2.0, 4.0]);
    /// assert_eq!(0.1 * vec3, Vector3::<f64>([0.2, 0.2, 0.4]));
    /// ```
    fn mul(self, rhs: Vector3<f64>) -> Self::Output {
        let v = rhs.map(|x| x * self);
        Vector3::<f64>(v)
    }
}

impl From<Vector3<Angstrom>> for Vector3<f64> {
    fn from(v: Vector3<Angstrom>) -> Self {
        Vector3::<f64>(v.map(f64::from))
    }
}

impl From<Vector3<f64>> for Vector3<Angstrom> {
    fn from(v: Vector3<f64>) -> Self {
        Vector3::<Angstrom>(v.map(Angstrom::from))
    }
}

impl From<Vector3<InvAngstrom>> for Vector3<f64> {
    fn from(v: Vector3<InvAngstrom>) -> Self {
        Vector3::<f64>(v.map(f64::from))
    }
}

impl From<Vector3<f64>> for Vector3<InvAngstrom> {
    fn from(v: Vector3<f64>) -> Self {
        Vector3::<InvAngstrom>(v.map(InvAngstrom::from))
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn matrix_3x3() {
        let mat = matrix_3x3![
            1 2 3;
            4 5 6.1;
            7 8 9;
        ];

        assert_eq!(mat[0][0], 1.0);
        assert_eq!(mat[1][2], 6.1);
    }
}
