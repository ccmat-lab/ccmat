use std::ops::{Deref, Mul};

use crate::structure::{Angstrom, InvAngstrom};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vec3<T>(pub [T; 3]);

impl<T> Deref for Vec3<T> {
    type Target = [T; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> std::ops::Index<usize> for Vec3<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Mul<Vec3<f64>> for f64 {
    type Output = Vec3<f64>;

    /// Scalar multiply for a `Vec3<f64>`
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmat_core::math::Vec3;
    ///
    /// let vec3 = Vec3::<f64>([2.0, 2.0, 4.0]);
    /// assert_eq!(0.1 * vec3, Vec3::<f64>([0.2, 0.2, 0.4]));
    /// ```
    fn mul(self, rhs: Vec3<f64>) -> Self::Output {
        let v = rhs.map(|x| x * self);
        Vec3::<f64>(v)
    }
}

impl From<Vec3<Angstrom>> for Vec3<f64> {
    fn from(v: Vec3<Angstrom>) -> Self {
        Vec3::<f64>(v.map(f64::from))
    }
}

impl From<Vec3<InvAngstrom>> for Vec3<f64> {
    fn from(v: Vec3<InvAngstrom>) -> Self {
        Vec3::<f64>(v.map(f64::from))
    }
}
