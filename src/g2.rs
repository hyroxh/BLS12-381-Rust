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

pub const orderG2: Fr = R_MODULUS;

pub const genG2: G2 = {
    G2 {
        valx: Fq2 {
            val0: Fq {
                limbs: [
                    0xc121bdb8, 0xd48056c8, 0xa805bbef, 0x0bac0326, 0x7ae3d177, 0xb4510b64,
                    0xfa403b02, 0xc6e47ad4, 0x2dc51051, 0x26080527, 0xf08f0a91, 0x024aa2b2,
                ],
            },

            val1: Fq {
                limbs: [
                    0x5d042b7e, 0xe5ac7d05, 0x13945d57, 0x334cf112, 0xdc7f5049, 0xb5da61bb,
                    0x9920b61a, 0x596bd0d0, 0x88274f65, 0x7dacd3a0, 0x52719f60, 0x13e02b60,
                ],
            },
        },

        valy: Fq2 {
            val0: Fq {
                limbs: [
                    0x08b82801, 0xe1935486, 0x3baca289, 0x923ac9cc, 0x5160d12c, 0x6d429a69,
                    0x8cbdd3a7, 0xadfd9baa, 0xda2e351a, 0x8cc9cdc6, 0x727d6e11, 0x0ce5d527,
                ],
            },

            val1: Fq {
                limbs: [
                    0xf05f79be, 0xaaa9075f, 0x5cec1da1, 0x3f370d27, 0x572e99ab, 0x267492ab,
                    0x85a763af, 0xcb3e287e, 0x2bc28b99, 0x32acd2b0, 0x2ea734cc, 0x0606c4a0,
                ],
            },
        },

        inf: 0u32,
    }
};

#[derive(Clone, Copy)]
pub struct G2 {
    pub valx: Fq2,
    pub valy: Fq2,
    pub inf: u32,
}

#[derive(Clone, Copy)]
pub struct G2Z {
    pub valX: Fq2,
    pub valY: Fq2,
    pub valZ: Fq2,
}

impl G2 {
    pub fn set(x: Fq2, y: Fq2) -> Self {
        Self {
            valx: x,
            valy: y,
            inf: 0u32,
        }
    }

    pub fn set_inf() -> Self {
        Self {
            valx: Fq2::zero(),
            valy: Fq2::zero(),
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

    pub fn to_project(&self) -> G2Z {
        //if self.inf = 0, mask == 0xF...F; if self.inf = 1, mask = 0
        let mask = !self.inf.wrapping_neg();

        let inf = G2Z::set_inf();
        let mut rop = G2Z::set_inf();

        for i in 0..N_LIMBS {
            rop.valX.val0.limbs[i] = self.valx.val0.limbs[i] & mask;
            rop.valX.val1.limbs[i] = self.valx.val1.limbs[i] & mask;

            rop.valY.val0.limbs[i] =
                (self.valy.val0.limbs[i] & mask) | (inf.valY.val0.limbs[i] & !mask);
            rop.valY.val1.limbs[i] =
                (self.valy.val1.limbs[i] & mask) | (inf.valY.val1.limbs[i] & !mask);

            rop.valZ.val0.limbs[i] = 0u32;
            rop.valZ.val1.limbs[i] = 0u32;
        }

        rop.valZ.val0.limbs[0] = (1u32 & mask) | (inf.valZ.val0.limbs[0] & !mask);
        rop
    }
}

impl G2Z {
    pub fn set(X: Fq2, Y: Fq2, Z: Fq2) -> Self {
        Self {
            valX: X,
            valY: Y,
            valZ: Z,
        }
    }

    pub fn set_inf() -> Self {
        Self {
            valX: Fq2::zero(),
            valY: Fq2::one(),
            valZ: Fq2::zero(),
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

    pub fn to_affine(&self) -> G2 {
        let zero = Fq2::zero();

        let mut mask = Fq2::is_equal(self.valZ, zero); //mask = 1, if Z = 0; 0, otherwise
        mask = mask.wrapping_neg(); // mask = 0xF...F, if Z = 0; 0, otherwise

        let Z_inv = self.valZ.mont_mul_inv();

        let res1 = G2::set_inf();
        let res2 = G2 {
            valx: self.valX * Z_inv,
            valy: self.valY * Z_inv,
            inf: 0u32,
        };

        let mut rop = G2::set_inf();

        for i in 0..N_LIMBS {
            rop.valx.val0.limbs[i] =
                (res1.valx.val0.limbs[i] & mask) | (res2.valx.val0.limbs[i] & !mask);
            rop.valx.val1.limbs[i] =
                (res1.valx.val1.limbs[i] & mask) | (res2.valx.val1.limbs[i] & !mask);

            rop.valy.val0.limbs[i] =
                (res1.valy.val0.limbs[i] & mask) | (res2.valy.val0.limbs[i] & !mask);
            rop.valy.val1.limbs[i] =
                (res1.valy.val1.limbs[i] & mask) | (res2.valy.val1.limbs[i] & !mask);
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

        mask &= Fq2::is_equal(XpYq, XqYp);
        mask &= Fq2::is_equal(YpZq, YqZp);
        mask &= Fq2::is_equal(ZpXq, ZqXp);

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
        let c8 = Fq2::new_ui(8u32, 0u32).mont_rep();
        let c18 = Fq2::new_ui(18u32, 0u32).mont_rep();
        let c16 = Fq2::new_ui(16u32, 0u32).mont_rep();
        let c27 = Fq2::new_ui(27u32, 0u32).mont_rep();
        let c36 = Fq2::new_ui(36u32, 0u32).mont_rep();

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
        let mut R0 = G2Z::set_inf();
        let mut R1 = self;
        let mut T_swap = G2Z::set_inf();

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
                T_swap.valX.val0.limbs[j] = (R0.valX.val0.limbs[j] ^ R1.valX.val0.limbs[j]) & mask;
                T_swap.valX.val1.limbs[j] = (R0.valX.val1.limbs[j] ^ R1.valX.val1.limbs[j]) & mask;
                T_swap.valY.val0.limbs[j] = (R0.valY.val0.limbs[j] ^ R1.valY.val0.limbs[j]) & mask;
                T_swap.valY.val1.limbs[j] = (R0.valY.val1.limbs[j] ^ R1.valY.val1.limbs[j]) & mask;
                T_swap.valZ.val0.limbs[j] = (R0.valZ.val0.limbs[j] ^ R1.valZ.val0.limbs[j]) & mask;
                T_swap.valZ.val1.limbs[j] = (R0.valZ.val1.limbs[j] ^ R1.valZ.val1.limbs[j]) & mask;

                R0.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
                R0.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
                R0.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
                R0.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
                R0.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
                R0.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];

                R1.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
                R1.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
                R1.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
                R1.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
                R1.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
                R1.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];
            }

            R1 = R0 + R1;
            R0 = R0 + R0;

            //swapping back
            for j in 0..N_LIMBS {
                T_swap.valX.val0.limbs[j] = (R0.valX.val0.limbs[j] ^ R1.valX.val0.limbs[j]) & mask;
                T_swap.valX.val1.limbs[j] = (R0.valX.val1.limbs[j] ^ R1.valX.val1.limbs[j]) & mask;
                T_swap.valY.val0.limbs[j] = (R0.valY.val0.limbs[j] ^ R1.valY.val0.limbs[j]) & mask;
                T_swap.valY.val1.limbs[j] = (R0.valY.val1.limbs[j] ^ R1.valY.val1.limbs[j]) & mask;
                T_swap.valZ.val0.limbs[j] = (R0.valZ.val0.limbs[j] ^ R1.valZ.val0.limbs[j]) & mask;
                T_swap.valZ.val1.limbs[j] = (R0.valZ.val1.limbs[j] ^ R1.valZ.val1.limbs[j]) & mask;

                R0.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
                R0.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
                R0.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
                R0.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
                R0.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
                R0.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];

                R1.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
                R1.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
                R1.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
                R1.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
                R1.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
                R1.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];
            }
        }

        R0
    }
}

// --- Arithmetic Trait Implementations ---

impl Add for G2Z {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let zero = Fq2::zero();

        let mask_inf1 = Fq2::is_equal(self.valZ, zero).wrapping_neg(); //0xF...F, if p=0; 0, otherwise
        let mask_inf2 = Fq2::is_equal(other.valZ, zero).wrapping_neg(); //0xF...F, if q=0; 0, otherwise
        let mask_inf = !mask_inf1 & !mask_inf2; //0xF...F, if p!=0 & q!=0; 0, otherwise

        let res_add = G2Z::adding(self, other);
        let res_dbl = G2Z::doubling(self);

        let mask = G2Z::is_equal(self, other).wrapping_neg(); //0xF..F, if p = q; 0, otherwise

        let mut rop = G2Z::set_inf();

        for i in 0..N_LIMBS {
            rop.valX.val0.limbs[i] = (((res_add.valX.val0.limbs[i] & !mask)
                | (res_dbl.valX.val0.limbs[i] & mask))
                & mask_inf)
                | (self.valX.val0.limbs[i] & mask_inf2)
                | (other.valX.val0.limbs[i] & mask_inf1);
            rop.valX.val1.limbs[i] = (((res_add.valX.val1.limbs[i] & !mask)
                | (res_dbl.valX.val1.limbs[i] & mask))
                & mask_inf)
                | (self.valX.val1.limbs[i] & mask_inf2)
                | (other.valX.val1.limbs[i] & mask_inf1);

            rop.valY.val0.limbs[i] = (((res_add.valY.val0.limbs[i] & !mask)
                | (res_dbl.valY.val0.limbs[i] & mask))
                & mask_inf)
                | (self.valY.val0.limbs[i] & mask_inf2)
                | (other.valY.val0.limbs[i] & mask_inf1);
            rop.valY.val1.limbs[i] = (((res_add.valY.val1.limbs[i] & !mask)
                | (res_dbl.valY.val1.limbs[i] & mask))
                & mask_inf)
                | (self.valY.val1.limbs[i] & mask_inf2)
                | (other.valY.val1.limbs[i] & mask_inf1);

            rop.valZ.val0.limbs[i] = (((res_add.valZ.val0.limbs[i] & !mask)
                | (res_dbl.valZ.val0.limbs[i] & mask))
                & mask_inf)
                | (self.valZ.val0.limbs[i] & mask_inf2)
                | (other.valZ.val0.limbs[i] & mask_inf1);
            rop.valZ.val1.limbs[i] = (((res_add.valZ.val1.limbs[i] & !mask)
                | (res_dbl.valZ.val1.limbs[i] & mask))
                & mask_inf)
                | (self.valZ.val1.limbs[i] & mask_inf2)
                | (other.valZ.val1.limbs[i] & mask_inf1);
        }

        rop
    }
}

impl Neg for G2Z {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            valX: self.valX,
            valY: -self.valY,
            valZ: self.valZ,
        }
    }
}

impl Sub for G2Z {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl Mul<Fr> for G2Z {
    type Output = Self;

    fn mul(self, other: Fr) -> Self {
        self.mul(other)
    }
}

impl Mul<G2Z> for Fr {
    type Output = G2Z;

    fn mul(self, other: G2Z) -> G2Z {
        other * self
    }
}

// --- Display / Printing ---

impl fmt::Display for G2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "<{}, {}, {}>", self.valx, self.valy, self.inf)
    }
}

impl fmt::Display for G2Z {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "<{} : {} : {}>", self.valX, self.valY, self.valZ)
    }
}
