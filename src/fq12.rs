use crate::fq::Fq;
use crate::fq2::Fq2;
use crate::fq6::Fq6;

use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy)]
pub struct Fq12 {
    pub val0: Fq6,
    pub val1: Fq6,
}

impl Fq12 {
    pub fn set(c0: Fq6, c1: Fq6) -> Self {
        Self { val0: c0, val1: c1 }
    }

    pub fn one() -> Self {
        Self {
            val0: Fq6::one(),
            val1: Fq6::zero(),
        }
    }

    pub fn zero() -> Self {
        Self {
            val0: Fq6::zero(),
            val1: Fq6::zero(),
        }
    }

    pub fn conj(&self) -> Self {
        Self {
            val0: self.val0,
            val1: -self.val1,
        }
    }
    pub fn is_equal(self, other: Self) -> u32 {
        Fq6::is_equal(self.val0, other.val0) & Fq6::is_equal(self.val1, other.val1)
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

    /*

    a*b = a0b0 + (a0b1 + a1b0)w + a1b1 w^2 = [w^2 = v] =

        = (a0b0 + a1b1 * v) + (a0b1 + a1b0)w.

    */

    fn mont_mul(&self, other: &Self) -> Self {
        let Self { val0: a0, val1: a1 } = *self;
        let Self { val0: b0, val1: b1 } = *other;

        let var00 = a0 * b0;
        let var01 = a0 * b1;
        let var10 = a1 * b0;
        let var11 = a1 * b1;

        let var11_t = var11.transform();

        Self {
            val0: var00 + var11_t,
            val1: var01 + var10,
        }
    }

    /*

    a^(-1) = (a0^2 - a1^2 * v)^(-1) * a0 + (a0^2 - a1^2 * v)^(-1) * (-a1) * w.

    */

    pub fn mont_mul_inv(&self) -> Self {
        let Self { val0: a0, val1: a1 } = *self;
        let var00 = a0 * a0;
        let var11 = a1 * a1;
        let var11_t = var11.transform();
        let norm = var00 - var11_t;
        let norm_inv = norm.mont_mul_inv();

        Self {
            val0: norm_inv * a0,
            val1: norm_inv * (-a1),
        }
    }
}

// --- Arithmetic Trait Implementations ---

impl Add for Fq12 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            val0: self.val0 + other.val0,
            val1: self.val1 + other.val1,
        }
    }
}

impl Sub for Fq12 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            val0: self.val0 - other.val0,
            val1: self.val1 - other.val1,
        }
    }
}

impl Neg for Fq12 {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            val0: -self.val0,
            val1: -self.val1,
        }
    }
}

impl Mul for Fq12 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.mont_mul(&other)
    }
}

// --- Display / Printing ---

impl fmt::Display for Fq12 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{{}, {}}}", self.val0, self.val1)
    }
}
