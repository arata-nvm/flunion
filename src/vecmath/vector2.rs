use std::{
    fmt,
    ops::{Deref, DerefMut},
};

#[derive(Debug, Default, Clone, Copy)]
pub struct Vector2(pub [f64; 2]);

impl Vector2 {
    pub fn new(x: f64, y: f64) -> Self {
        Self([x, y])
    }

    pub fn zero() -> Self {
        Self::default()
    }

    pub fn set_zero(&mut self) {
        self[0] = 0.0;
        self[1] = 0.0;
    }

    pub fn set(&mut self, x: f64, y: f64) {
        self[0] = x;
        self[1] = y;
    }

    pub fn add(&self, v: Vector2) -> Self {
        Self::new(self[0] + v[0], self[1] + v[1])
    }

    pub fn sub(&self, v: Vector2) -> Self {
        Self::new(self[0] - v[0], self[1] - v[1])
    }

    pub fn mul(&self, f: f64) -> Self {
        Self::new(self[0] * f, self[1] * f)
    }

    pub fn div(&self, f: f64) -> Self {
        if f == 0.0 {
            return Self::zero();
        }

        Self::new(self[0] / f, self[1] / f)
    }

    pub fn set_add(&mut self, v: Vector2) {
        self[0] += v[0];
        self[1] += v[1];
    }

    pub fn set_sub(&mut self, v: Vector2) {
        self[0] -= v[0];
        self[1] -= v[1];
    }

    pub fn set_mul(&mut self, f: f64) {
        self[0] *= f;
        self[1] *= f;
    }

    pub fn set_div(&mut self, f: f64) {
        if f == 0.0 {
            self.set_zero();
        } else {
            self[0] /= f;
            self[1] /= f;
        }
    }

    pub fn len(&self) -> f64 {
        f64::sqrt(self[0] * self[0] + self[1] * self[1])
    }

    pub fn len_sq(&self) -> f64 {
        self[0] * self[0] + self[1] * self[1]
    }

    pub fn normalized(&self) -> Self {
        self.div(self.len())
    }

    pub fn neg(&self) -> Self {
        Self::new(-self[0], -self[1])
    }

    pub fn rot(&self, r: f64) -> Self {
        let c = f64::cos(r);
        let s = f64::sin(r);
        Self::new(c * self[0] - s * self[1], s * self[0] + c * self[1])
    }

    pub fn rot_rev(&self, r: f64) -> Self {
        let c = f64::cos(r);
        let s = f64::sin(r);
        Self::new(c * self[0] + s * self[1], -s * self[0] + c * self[1])
    }
}

pub fn dot(v1: Vector2, v2: Vector2) -> f64 {
    v1[0] * v2[0] + v1[1] * v2[1]
}

pub fn cross(v1: Vector2, v2: Vector2) -> f64 {
    v1[0] * v2[1] - v1[1] * v2[0]
}

pub fn crossf(v: Vector2, f: f64) -> Vector2 {
    Vector2::new(-f * v[1], v[0] * f)
}

impl Deref for Vector2 {
    type Target = [f64; 2];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Vector2 {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Vector2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self[0], self[1])
    }
}
