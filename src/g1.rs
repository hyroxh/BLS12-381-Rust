use crate::fq::Fq;
use crate::fq2::Fq2;
use crate::fq6::Fq6;
use crate::fq12::Fq12;
use crate::fr::Fr;

use crate::fq::N_LIMBS;

use crate::fr::FR_LIMBS;
use crate::fr::LEN_MAX_R;
use crate::fr::R_MODULUS;

use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

pub const orderG1: Fr = R_MODULUS;

pub const genG1: G1 = {
    G1 {
        valx: Fq {
            limbs: [
                0xdb22c6bb, 0xfb3af00a, 0xf97a1aef, 0x6c55e83f, 0x171bac58, 0xa14e3a3f, 0x9774b905,
                0xc3688c4f, 0x4fa9ac0f, 0x2695638c, 0x3197d794, 0x17f1d3a7,
            ],
        },

        valy: Fq {
            limbs: [
                0x46c5e7e1, 0x0caa2329, 0xa2888ae4, 0xd03cc744, 0x2c04b3ed, 0x00db18cb, 0xd5d00af6,
                0xfcf5e095, 0x741d8ae4, 0xa09e30ed, 0xe3aaa0f1, 0x08b3f481,
            ],
        },

        inf: 0u32,
    }
};

#[derive(Clone, Copy)]
pub struct G1 {
    pub valx: Fq,
    pub valy: Fq,
    pub inf: u32,
}

#[derive(Clone, Copy)]
pub struct G1Z {
    pub valX: Fq,
    pub valY: Fq,
    pub valZ: Fq,
}

impl G1 {
    pub fn set(x: Fq, y: Fq) -> Self {
        Self {
            valx: x,
            valy: y,
            inf: 0u32,
        }
    }

    pub fn set_inf() -> Self {
        Self {
            valx: Fq::zero(),
            valy: Fq::zero(),
            inf: 1u32,
        }
    }

    pub fn mont_rep(&self) -> Self {
        Self {
            valx: self.valx.mont_rep(),
            valy: self.valy.mont_rep(),
            inf: self.inf,
        }
    }

    pub fn mont_rep_inv(&self) -> Self {
        Self {
            valx: self.valx.mont_rep_inv(),
            valy: self.valy.mont_rep_inv(),
            inf: self.inf,
        }
    }

    pub fn to_project(&self) -> G1Z {
        //if self.inf = 0, mask == 0xF...F; if self.inf = 1, mask = 0
        let mask = !self.inf.wrapping_neg();

        let inf = G1Z::set_inf();
        let mut rop = G1Z::set_inf();

        for i in 0..N_LIMBS {
            rop.valX.limbs[i] = self.valx.limbs[i] & mask;
            rop.valY.limbs[i] = (self.valy.limbs[i] & mask) | (inf.valY.limbs[i] & !mask);
            rop.valZ.limbs[i] = 0u32;
        }

        rop.valZ.limbs[0] = (1u32 & mask) | (inf.valZ.limbs[0] & !mask);
        rop
    }
}

impl G1Z {
    pub fn set(X: Fq, Y: Fq, Z: Fq) -> Self {
        Self {
            valX: X,
            valY: Y,
            valZ: Z,
        }
    }

    pub fn set_inf() -> Self {
        Self {
            valX: Fq::zero(),
            valY: Fq::one(),
            valZ: Fq::zero(),
        }
    }

    pub fn mont_rep(&self) -> Self {
        Self {
            valX: self.valX.mont_rep(),
            valY: self.valY.mont_rep(),
            valZ: self.valZ.mont_rep(),
        }
    }

    pub fn mont_rep_inv(&self) -> Self {
        Self {
            valX: self.valX.mont_rep_inv(),
            valY: self.valY.mont_rep_inv(),
            valZ: self.valZ.mont_rep_inv(),
        }
    }

    pub fn to_affine(&self) -> G1 {
        let zero = Fq::zero();

        let mut mask = Fq::is_equal(self.valZ, zero); //mask = 1, if Z = 0; 0, otherwise
        mask = mask.wrapping_neg(); // mask = 0xF...F, if Z = 0; 0, otherwise

        let Z_inv = self.valZ.mont_mul_inv();

        let res1 = G1::set_inf();
        let res2 = G1 {
            valx: self.valX * Z_inv,
            valy: self.valY * Z_inv,
            inf: 0u32,
        };

        let mut rop = G1::set_inf();

        for i in 0..N_LIMBS {
            rop.valx.limbs[i] = (res1.valx.limbs[i] & mask) | (res2.valx.limbs[i] & !mask);
            rop.valy.limbs[i] = (res1.valy.limbs[i] & mask) | (res2.valy.limbs[i] & !mask);
        }

        rop.inf = (res1.inf & mask) | (res2.inf & !mask);
        rop
    }

    //return 1, if P = Q; 0, otherwise
    pub fn is_equal(self, other: Self) -> u32 {
        //(Xp : Yp : Zp) ~ (Xq : Yq : Zq) <=>
        // <=> [XpYq - XqYp, YpZq - YqZp, ZpXq - ZqXp] = [0, 0, 0]

        let mut mask = 1u32;

        let XpYq = self.valX * other.valY;
        let XqYp = self.valY * other.valX;
        let YpZq = self.valY * other.valZ;
        let YqZp = self.valZ * other.valY;
        let ZpXq = self.valZ * other.valX;
        let ZqXp = self.valX * other.valZ;

        mask &= Fq::is_equal(XpYq, XqYp);
        mask &= Fq::is_equal(YpZq, YqZp);
        mask &= Fq::is_equal(ZpXq, ZqXp);

        mask
    }

    // If p!=q
    fn adding(self, other: Self) -> Self {
        /*

        Zr = ZpZq * (XpZq - XqZp)^3.

        Xr = (XpZq - XqZp) * [ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq +
        XqZp)]

        Yr = ZpZq * (XqYp - XpYq) * (XpZq - XqZp)^2 - (YpZq - YqZp) * [ZpZq *
        (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

        */

        let Xp = self.valX;
        let Yp = self.valY;
        let Zp = self.valZ;

        let Xq = other.valX;
        let Yq = other.valY;
        let Zq = other.valZ;

        let Zr = Zp * Zq * (Xp * Zq - Xq * Zp) * (Xp * Zq - Xq * Zp) * (Xp * Zq - Xq * Zp);
        let Xr = (Xp * Zq - Xq * Zp)
            * (Zp * Zq * (Yp * Zq - Yq * Zp) * (Yp * Zq - Yq * Zp)
                - (Xp * Zq - Xq * Zp) * (Xp * Zq - Xq * Zp) * (Xp * Zq + Xq * Zp));
        let Yr = Zp * Zq * (Xq * Yp - Xp * Yq) * (Xp * Zq - Xq * Zp) * (Xp * Zq - Xq * Zp)
            - (Yp * Zq - Yq * Zp)
                * (Zp * Zq * (Yp * Zq - Yq * Zp) * (Yp * Zq - Yq * Zp)
                    - (Xp * Zq - Xq * Zp) * (Xp * Zq - Xq * Zp) * (Xp * Zq + Xq * Zp));

        Self {
            valX: Xr,
            valY: Yr,
            valZ: Zr,
        }
    }

    // If p=q
    fn doubling(self) -> Self {
        /*

        Zr = 8Yp^3Zp^3

        Xr = 18Xp^4YpZp - 16XpYp^3Zp^2

        Yr = -27Xp^6 + 36Xp^3Yp^2Zp - 8Yp^4Zp^2

        */
        let c8 = Fq::new_ui(8u32).mont_rep();
        let c18 = Fq::new_ui(18u32).mont_rep();
        let c16 = Fq::new_ui(16u32).mont_rep();
        let c27 = Fq::new_ui(27u32).mont_rep();
        let c36 = Fq::new_ui(36u32).mont_rep();

        let Xp = self.valX;
        let Yp = self.valY;
        let Zp = self.valZ;

        let Zr = c8 * Yp * Yp * Yp * Zp * Zp * Zp;
        let Xr = c18 * Xp * Xp * Xp * Xp * Yp * Zp - c16 * Xp * Yp * Yp * Yp * Zp * Zp;
        let Yr = c36 * Xp * Xp * Xp * Yp * Yp * Zp
            - c27 * Xp * Xp * Xp * Xp * Xp * Xp
            - c8 * Yp * Yp * Yp * Yp * Zp * Zp;

        Self {
            valX: Xr,
            valY: Yr,
            valZ: Zr,
        }
    }

    // rop = [n]P, Montgomery ladder
    fn mul(self, other: Fr) -> Self {
        let mut R0 = G1Z::set_inf();
        let mut R1 = self;
        let mut T_swap = G1Z::set_inf();

        let mut binstr = [0u32; 4 * LEN_MAX_R];

        for i in 0..FR_LIMBS {
            for j in 0..32 {
                binstr[32 * i + j] = (other.limbs[i] >> j) & 1u32;
            }
        }

        for i in (0..=4 * LEN_MAX_R - 1).rev() {
            let mask = binstr[i].wrapping_neg(); ////0xFFFFFFFF, if 1; 0, otherwise

            //swapping
            for j in 0..N_LIMBS {
                T_swap.valX.limbs[j] = (R0.valX.limbs[j] ^ R1.valX.limbs[j]) & mask;
                T_swap.valY.limbs[j] = (R0.valY.limbs[j] ^ R1.valY.limbs[j]) & mask;
                T_swap.valZ.limbs[j] = (R0.valZ.limbs[j] ^ R1.valZ.limbs[j]) & mask;

                R0.valX.limbs[j] ^= T_swap.valX.limbs[j];
                R0.valY.limbs[j] ^= T_swap.valY.limbs[j];
                R0.valZ.limbs[j] ^= T_swap.valZ.limbs[j];

                R1.valX.limbs[j] ^= T_swap.valX.limbs[j];
                R1.valY.limbs[j] ^= T_swap.valY.limbs[j];
                R1.valZ.limbs[j] ^= T_swap.valZ.limbs[j];
            }

            R1 = R0 + R1;
            R0 = R0 + R0;

            //swapping back
            for j in 0..N_LIMBS {
                T_swap.valX.limbs[j] = (R0.valX.limbs[j] ^ R1.valX.limbs[j]) & mask;
                T_swap.valY.limbs[j] = (R0.valY.limbs[j] ^ R1.valY.limbs[j]) & mask;
                T_swap.valZ.limbs[j] = (R0.valZ.limbs[j] ^ R1.valZ.limbs[j]) & mask;

                R0.valX.limbs[j] ^= T_swap.valX.limbs[j];
                R0.valY.limbs[j] ^= T_swap.valY.limbs[j];
                R0.valZ.limbs[j] ^= T_swap.valZ.limbs[j];

                R1.valX.limbs[j] ^= T_swap.valX.limbs[j];
                R1.valY.limbs[j] ^= T_swap.valY.limbs[j];
                R1.valZ.limbs[j] ^= T_swap.valZ.limbs[j];
            }
        }

        R0
    }
}

// --- Arithmetic Trait Implementations ---

impl Add for G1Z {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let zero = Fq::zero();

        let mask_inf1 = Fq::is_equal(self.valZ, zero).wrapping_neg(); //0xF...F, if p=0; 0, otherwise
        let mask_inf2 = Fq::is_equal(other.valZ, zero).wrapping_neg(); //0xF...F, if q=0; 0, otherwise
        let mask_inf = !mask_inf1 & !mask_inf2; //0xF...F, if p!=0 & q!=0; 0, otherwise

        let res_add = G1Z::adding(self, other);
        let res_dbl = G1Z::doubling(self);

        let mask = G1Z::is_equal(self, other).wrapping_neg(); //0xF..F, if p = q; 0, otherwise

        let mut rop = G1Z::set_inf();

        for i in 0..N_LIMBS {
            rop.valX.limbs[i] =
                (((res_add.valX.limbs[i] & !mask) | (res_dbl.valX.limbs[i] & mask)) & mask_inf)
                    | (self.valX.limbs[i] & mask_inf2)
                    | (other.valX.limbs[i] & mask_inf1);
            rop.valY.limbs[i] =
                (((res_add.valY.limbs[i] & !mask) | (res_dbl.valY.limbs[i] & mask)) & mask_inf)
                    | (self.valY.limbs[i] & mask_inf2)
                    | (other.valY.limbs[i] & mask_inf1);
            rop.valZ.limbs[i] =
                (((res_add.valZ.limbs[i] & !mask) | (res_dbl.valZ.limbs[i] & mask)) & mask_inf)
                    | (self.valZ.limbs[i] & mask_inf2)
                    | (other.valZ.limbs[i] & mask_inf1);
        }

        rop
    }
}

impl Neg for G1Z {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            valX: self.valX,
            valY: -self.valY,
            valZ: self.valZ,
        }
    }
}

impl Sub for G1Z {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl Mul<Fr> for G1Z {
    type Output = Self;

    fn mul(self, other: Fr) -> Self {
        self.mul(other)
    }
}

impl Mul<G1Z> for Fr {
    type Output = G1Z;

    fn mul(self, other: G1Z) -> G1Z {
        other * self
    }
}

// --- Display / Printing ---

impl fmt::Display for G1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "<{}, {}, {}>", self.valx, self.valy, self.inf)
    }
}

impl fmt::Display for G1Z {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "<{} : {} : {}>", self.valX, self.valY, self.valZ)
    }
}
