use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};
use std::ptr;

// Number of 32-bit limbs (12 * 32 = 384 bits)
pub const N_LIMBS: usize = 12;
pub const LEN_MAX: usize = 96;

// The modulus q = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
pub const Q_MODULUS: Fq = Fq {
    limbs: [
        0xffffaaab, 0xb9feffff, 0xb153ffff, 0x1eabfffe, 0xf6b0f624, 0x6730d2a0, 0xf38512bf,
        0x64774b84, 0x434bacd7, 0x4b1ba7b6, 0x397fe69a, 0x1a0111ea,
    ],
};

// Q_LSLIMB_INV = -(q^{-1} mod 2^32)
pub const Q_LSLIMB_INV: u32 = 196611;

// Binary representation of (q - 2) for inversion
// Used in Fq_mont_mul_inv
pub const Q_MINUS_2_BIN: &str = "110100000000100010001111010100011100101111111111\
0011010011010010010110001101110100111101101100100001101001011\
10101100110101110110010001110111010010111000010011110011100001\
01000100101011111101100111001100001101001010100000111101101011\
00001111011000100100000111101010101111111111111111101011000101\
01001111111111111111111011100111111110111111111111111111111111\
111111111010101010101001";

#[derive(Clone, Copy)]
pub struct Fq {
    pub limbs: [u32; N_LIMBS],
}

impl Fq {
    // Create a new Fq from a u32 (sets LS limb)
    pub fn new_ui(val: u32) -> Self {
        let mut limbs = [0u32; N_LIMBS];
        limbs[0] = val;
        Fq { limbs }
    }

    pub fn zero() -> Self {
        Fq {
            limbs: [0u32; N_LIMBS],
        }
    }

    pub fn one() -> Self {
        Self::new_ui(1)
    }

    //return 1, if a = b; 0, otherwise
    pub fn is_equal(self, other: Self) -> u32 {
        let mut res = 0u32;
        for i in 0..N_LIMBS {
            res |= self.limbs[i] ^ other.limbs[i];
        }

        //if a = b, res = 0; if a!=b, res = 1
        res = (res | res.wrapping_neg()) >> 31;

        //if a = b, res = 1; if a!=b, res = 0
        res = res ^ 1u32;

        res
    }

    pub fn from_hex_str(hex_str: &str) -> Result<Self, &'static str> {
        let mut limbs = [0u32; N_LIMBS];
        let clean_str = hex_str.trim_start_matches("0x");

        // We process the string from right to left (LSB to MSB)
        let chars: Vec<char> = clean_str.chars().collect();
        let len = chars.len();

        for i in 0..N_LIMBS {
            let mut val = 0u32;
            // Each limb is 8 hex digits
            for j in 0..8 {
                let char_idx = (i * 8) + j;
                if char_idx < len {
                    let c = chars[len - 1 - char_idx];
                    let nibble = c.to_digit(16).ok_or("Invalid hex character")?;
                    val |= nibble << (j * 4);
                }
            }
            limbs[i] = val;
        }

        Ok(Fq { limbs })
    }
    // If a == q, make a = 0; otherwise, leave untouched.
    // Constant-time implementation.
    fn correction(&mut self) {
        let mut mask = 0u32;
        for i in 0..N_LIMBS {
            mask |= self.limbs[i] ^ Q_MODULUS.limbs[i];
        }

        //if a == q, mask = 0; if a!=q, mask = 1
        mask = (mask | mask.wrapping_neg()) >> 31;
        // if a == q, mask = 0; if a!=q, mask = 1...1
        mask = mask.wrapping_neg();

        for i in 0..N_LIMBS {
            self.limbs[i] &= mask;
        }
    }

    // Montgomery Representation: aR mod q
    pub fn mont_rep(&self) -> Self {
        let mut aux = *self;
        // R = 2^(32 * N_LIMBS). We add R times (doubling)
        for _ in 0..(32 * N_LIMBS) {
            aux = aux + aux;
        }
        aux.correction();
        aux
    }

    // Convert back from Montgomery: aR^-1 mod q
    pub fn mont_rep_inv(&self) -> Self {
        let id = Fq::one();
        // mont_mul calculates (a * b * R^-1).
        // passing '1' results in (a * 1 * R^-1) = aR^-1
        self.mont_mul(&id)
    }

    // Montgomery Multiplication
    // Returns: a * b * R^-1 mod q
    fn mont_mul(&self, other: &Self) -> Self {
        // g = -Q_inv mod 2^32
        let g = (Q_LSLIMB_INV as i32).wrapping_neg() as u32;
        let mut d = Fq::zero();
        let mut dh: u32 = 0;

        for i in 0..N_LIMBS {
            // f = (d[0] + a[i]*b[0]) * g
            let a_limb = self.limbs[i];
            let term = d.limbs[0].wrapping_add(a_limb.wrapping_mul(other.limbs[0]));
            let f = term.wrapping_mul(g);

            let mut c: u64 = 0;

            for j in 0..N_LIMBS {
                // z = d[j] + a[i]*b[j] + f*q[j] + c
                let term1 = d.limbs[j] as u128;
                let term2 = (a_limb as u128) * (other.limbs[j] as u128);
                let term3 = (f as u128) * (Q_MODULUS.limbs[j] as u128);
                let z = term1 + term2 + term3 + (c as u128);

                if j > 0 {
                    d.limbs[j - 1] = z as u32;
                }
                c = (z >> 32) as u64;
            }

            let z = (dh as u64) + c;
            d.limbs[N_LIMBS - 1] = z as u32;
            dh = (z >> 32) as u32;
        }

        // Subtraction: d_sub_Q = d - Q
        let mut d_sub_q = Fq::zero();
        let mut borrow: u64 = 0;
        for i in 0..N_LIMBS {
            let diff = (d.limbs[i] as u64)
                .wrapping_sub(Q_MODULUS.limbs[i] as u64)
                .wrapping_sub(borrow);
            d_sub_q.limbs[i] = diff as u32;
            borrow = (diff >> 63) & 1;
        }

        //if dh >= borrow, mask_cond is 1; 0, otherwise
        let mask_cond = (((dh as u64).wrapping_sub(borrow as u64) >> 63) as u32) ^ 1;
        let mask = (mask_cond as i32).wrapping_neg() as u32;

        let mut rop = Fq::zero();
        for i in 0..N_LIMBS {
            rop.limbs[i] = (d_sub_q.limbs[i] & mask) | (d.limbs[i] & !mask);
        }

        rop.correction();
        rop
    }

    // Modular Inverse using Fermat's Little Theorem
    pub fn mont_mul_inv(&self) -> Self {
        let mut aux = *self;
        let chars: Vec<char> = Q_MINUS_2_BIN.chars().collect();

        for i in 1..chars.len() {
            aux = aux.mont_mul(&aux);
            if chars[i] == '1' {
                aux = aux.mont_mul(self);
            }
        }
        aux
    }
}

// --- Arithmetic Trait Implementations ---

impl Add for Fq {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut d = Fq::zero();
        let mut e = Fq::zero();
        let mut c1: u64 = 0;
        let mut c2: u64 = 0;

        // d = a + b
        for i in 0..N_LIMBS {
            let sum = (self.limbs[i] as u64) + (other.limbs[i] as u64) + c1;
            d.limbs[i] = sum as u32;
            c1 = sum >> 32;
        }

        // e = d - Q
        for i in 0..N_LIMBS {
            let diff = (d.limbs[i] as u64)
                .wrapping_sub(Q_MODULUS.limbs[i] as u64)
                .wrapping_sub(c2);
            e.limbs[i] = diff as u32;
            c2 = (diff >> 63) & 1;
        }

        // If c2 == 0, it means d >= Q, so we choose e.
        // If c2 == 1, it means d < Q, so we choose d.
        // mask = -(c2 == 0) -> 0xFFFFFFFF if c2=0
        let mask_cond = c2 ^ 1u64;
        let mask = (mask_cond as u32).wrapping_neg();

        let mut rop = Fq::zero();
        for i in 0..N_LIMBS {
            rop.limbs[i] = (e.limbs[i] & mask) | (d.limbs[i] & !mask);
        }
        rop.correction();
        rop
    }
}

impl Sub for Fq {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut d = Fq::zero();
        let mut e = Fq::zero();
        let mut c1: u64 = 0;
        let mut c2: u64 = 0;

        // d = a - b
        for i in 0..N_LIMBS {
            let diff = (self.limbs[i] as u64)
                .wrapping_sub(other.limbs[i] as u64)
                .wrapping_sub(c1);
            d.limbs[i] = diff as u32;
            c1 = (diff >> 63) & 1;
        }

        // e = d + Q
        for i in 0..N_LIMBS {
            let sum = (d.limbs[i] as u64) + (Q_MODULUS.limbs[i] as u64) + c2;
            e.limbs[i] = sum as u32;
            c2 = sum >> 32;
        }

        // if c1 == 1 (borrow occurred), result is negative, so we add Q (select e).
        // if c1 == 0, result is positive, we select d.
        let mask_cond = c1 as u32;
        let mask = mask_cond.wrapping_neg();

        let mut rop = Fq::zero();
        for i in 0..N_LIMBS {
            rop.limbs[i] = (e.limbs[i] & mask) | (d.limbs[i] & !mask);
        }
        rop.correction();
        rop
    }
}

impl Neg for Fq {
    type Output = Self;
    fn neg(self) -> Self {
        let q = Q_MODULUS;
        q - self
    }
}

impl Mul for Fq {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self.mont_mul(&other)
    }
}

// --- Display / Printing ---

impl fmt::Display for Fq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "0x")?;
        // Print from Most Significant Limb to Least
        for i in (0..N_LIMBS).rev() {
            write!(f, "{:08x}", self.limbs[i])?;
        }
        Ok(())
    }
}
