use crate::fq::Fq;

use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy)]
pub struct Fq2 {
    pub val0: Fq,
    pub val1: Fq,
}

impl Fq2 {
    pub fn new_ui(c0: u32, c1: u32) -> Self {
        Self {
            val0: Fq::new_ui(c0),
            val1: Fq::new_ui(c1),
        }
    }

    pub fn set(c0: Fq, c1: Fq) -> Self {
        Self { val0: c0, val1: c1 }
    }

    pub fn one() -> Self {
        Self {
            val0: Fq::one(),
            val1: Fq::zero(),
        }
    }

    pub fn zero() -> Self {
        Self {
            val0: Fq::zero(),
            val1: Fq::zero(),
        }
    }

    pub fn conj(&self) -> Self {
        Self {
            val0: self.val0,
            val1: -self.val1,
        }
    }

    pub fn is_equal(self, other: Self) -> u32 {
        Fq::is_equal(self.val0, other.val0) & Fq::is_equal(self.val1, other.val1)
    }

    pub fn mont_rep(&self) -> Self {
        Self {
            val0: self.val0.mont_rep(),
            val1: self.val1.mont_rep(),
        }
    }

    pub fn mont_rep_inv(&self) -> Self {
        Self {
            val0: self.val0.mont_rep_inv(),
            val1: self.val1.mont_rep_inv(),
        }
    }

    pub fn transform(&self) -> Self {
        Self {
            val0: -self.val1,
            val1: self.val0,
        }
    }

    fn mont_mul(&self, other: &Self) -> Self {
        let a = self.val0 * other.val0 - self.val1 * other.val1;
        let b = self.val0 * other.val1 + self.val1 * other.val0;

        Self { val0: a, val1: b }
    }

    pub fn mont_mul_inv(&self) -> Self {
        let norm = self.val0 * self.val0 + self.val1 * self.val1;
        let norm_inv = Fq::mont_mul_inv(&norm);

        Self {
            val0: norm_inv * self.val0,
            val1: norm_inv * (-self.val1),
        }
    }
}

// --- Arithmetic Trait Implementations ---

impl Add for Fq2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            val0: self.val0 + other.val0,
            val1: self.val1 + other.val1,
        }
    }
}

impl Sub for Fq2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            val0: self.val0 - other.val0,
            val1: self.val1 - other.val1,
        }
    }
}

impl Neg for Fq2 {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            val0: -self.val0,
            val1: -self.val1,
        }
    }
}

impl Mul for Fq2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.mont_mul(&other)
    }
}

// --- Display / Printing ---

impl fmt::Display for Fq2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.val0, self.val1)
    }
}
