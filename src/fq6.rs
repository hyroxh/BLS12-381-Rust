use crate::fq::Fq;
use crate::fq2::Fq2;

use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy)]
pub struct Fq6 {
    pub val0: Fq2,
    pub val1: Fq2,
    pub val2: Fq2,
}

impl Fq6 {
    pub fn set(c0: Fq2, c1: Fq2, c2: Fq2) -> Self {
        Self {
            val0: c0,
            val1: c1,
            val2: c2,
        }
    }

    pub fn one() -> Self {
        Self {
            val0: Fq2::one(),
            val1: Fq2::zero(),
            val2: Fq2::zero(),
        }
    }

    pub fn zero() -> Self {
        Self {
            val0: Fq2::zero(),
            val1: Fq2::zero(),
            val2: Fq2::zero(),
        }
    }
	
	pub fn is_equal(self, other: Self) -> u32 {
        Fq2::is_equal(self.val0, other.val0) & Fq2::is_equal(self.val1, other.val1) & Fq2::is_equal(self.val2, other.val2)
    }
    pub fn mont_rep(&self) -> Self {
        Self {
            val0: self.val0.mont_rep(),
            val1: self.val1.mont_rep(),
            val2: self.val2.mont_rep(),
        }
    }

    pub fn mont_rep_inv(&self) -> Self {
        Self {
            val0: self.val0.mont_rep_inv(),
            val1: self.val1.mont_rep_inv(),
            val2: self.val2.mont_rep_inv(),
        }
    }

    // a -> a * v
    pub fn transform(&self) -> Self {
        let u1 = Fq2::new_ui(1, 1);
        let u1m = u1.mont_rep();

        Self {
            val0: self.val2 * u1m,
            val1: self.val0,
            val2: self.val1,
        }
    }

    fn determinant(&self) -> Fq2 {
        /*

        det M = a0^3 - 3a0a1a2(u+1) + a1^3(u+1) + 2a2^3 u = Norm(a) = a * a^(q^2) *
            * a^(q^4).

        */
        let Self {
            val0: a0,
            val1: a1,
            val2: a2,
        } = *self;

        let u = Fq2::new_ui(0, 1);
        let um = u.mont_rep();

        let u1 = Fq2::new_ui(1, 1);
        let u1m = u1.mont_rep();

        let det_m = a0 * a0 * a0 - a0 * a1 * a2 * u1m - a0 * a1 * a2 * u1m - a0 * a1 * a2 * u1m
            + a1 * a1 * a1 * u1m
            + a2 * a2 * a2 * um
            + a2 * a2 * a2 * um;

        det_m
    }

    /*

    a*b = a0b0 + (a0b1 + a1b0)v + (a0b2 + a1b1 + a2b0)v^2 + (a1b2 + a2b1)v^3 +
        a2b2v^4 = [v^3 = u + 1] =

    = [a0b0 + (a1b2 + a2b1)(u + 1)] + [a0b1 + a1b0 + (a2b2)(u + 1)]v +
        [a0b2 + a1b1 + a2b0]v^2.

    */

    fn mont_mul(&self, other: &Self) -> Self {
        let u1 = Fq2::new_ui(1, 1);
        let u1m = u1.mont_rep();

        let Self {
            val0: a0,
            val1: a1,
            val2: a2,
        } = *self;
        let Self {
            val0: b0,
            val1: b1,
            val2: b2,
        } = *other;

        Self {
            val0: a0 * b0 + (a1 * b2 + a2 * b1) * u1m,
            val1: a0 * b1 + a1 * b0 + a2 * b2 * u1m,
            val2: a0 * b2 + a1 * b1 + a2 * b0,
        }
    }

    /*

    a^(-1) = detM ^(-1) * ((a0^2 - a1a2(u+1)) + (a2^2(u+1) - a0a1)v + (a1^2 -
        - a0a2)v^2).

    */

    pub fn mont_mul_inv(&self) -> Self {
        let det = self.determinant();
        let detinv = det.mont_mul_inv();
        let u1 = Fq2::new_ui(1, 1);
        let u1m = u1.mont_rep();

        let Self {
            val0: a0,
            val1: a1,
            val2: a2,
        } = *self;

        Self {
            val0: (a0 * a0 - a1 * a2 * u1m) * detinv,
            val1: (a2 * a2 * u1m - a0 * a1) * detinv,
            val2: (a1 * a1 - a0 * a2) * detinv,
        }
    }
}

// --- Arithmetic Trait Implementations ---

impl Add for Fq6 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            val0: self.val0 + other.val0,
            val1: self.val1 + other.val1,
            val2: self.val2 + other.val2,
        }
    }
}

impl Sub for Fq6 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            val0: self.val0 - other.val0,
            val1: self.val1 - other.val1,
            val2: self.val2 - other.val2,
        }
    }
}

impl Neg for Fq6 {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            val0: -self.val0,
            val1: -self.val1,
            val2: -self.val2,
        }
    }
}

impl Mul for Fq6 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.mont_mul(&other)
    }
}

// --- Display / Printing ---

impl fmt::Display for Fq6 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}, {}]", self.val0, self.val1, self.val2)
    }
}
