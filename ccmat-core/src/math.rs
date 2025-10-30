use std::ops::{Deref, Mul};

use crate::structure::{Angstrom, InvAngstrom};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vector3<T>(pub [T; 3]);

pub(crate) type TransformationMatrix = [[i32; 3]; 3];

impl<T> Deref for Vector3<T> {
    type Target = [T; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> std::ops::Index<usize> for Vector3<T> {
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
