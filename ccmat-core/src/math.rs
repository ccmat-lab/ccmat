use std::ops::{Deref, Mul};

use crate::structure::Angstrom;

#[derive(Debug, Copy, Clone)]
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

    fn mul(self, rhs: Vec3<f64>) -> Self::Output {
        let v = rhs.0.map(f64::from);
        Vec3::<f64>(v)
    }
}

impl From<Vec3<Angstrom>> for Vec3<f64> {
    fn from(v: Vec3<Angstrom>) -> Self {
        Vec3::<f64>(v.map(f64::from))
    }
}
