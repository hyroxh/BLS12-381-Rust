use crate::fq::Fq;
use crate::fq2::Fq2;
use crate::fq6::Fq6;
use crate::fq12::Fq12;
use crate::fr::Fr;
use crate::g1::G1;
use crate::g1::G1Z;
use crate::g2::G2;
use crate::g2::G2Z;

use crate::fq::N_LIMBS;

use crate::fr::FR_LIMBS;
use crate::fr::LEN_MAX_R;
use crate::fr::R_MODULUS;

use crate::g1::genG1;
use crate::g1::orderG1;
use crate::g2::genG2;
use crate::g2::orderG2;

use crate::frobenius;

const bATE_STR: &str = "1101001000000001000000000000000000000000000000010000000000000000"; //T = |x|, with x = t - 1, where t is the Frobenius trace
const len_bATE_STR: usize = 64;

fn ate_line_doubling(RZ: G2Z, P: G1) -> Fq12 {
    /*
    Computing L_R_R (P)

    L_R_R (P) = (2*YrZr^2) * Yp - 3*Xr^2Zr * Xp + 3Xr^3 - 2Yr^2Zr ->

    -> (2*YrZr^2) * Yp * const * vw - 3*Xr^2Zr * Xp * const * v + 3Xr^3 * const
    - 2Yr^2Zr * const,

    where const = sextic_const = (1-u)/2 <--- no need to multiply due to Frobenius endomorphism

    */

    let Xr = RZ.valX;
    let Yr = RZ.valY;
    let Zr = RZ.valZ;

    let zero = Fq::zero();

    let Xp = Fq2::set(P.valx, zero);
    let Yp = Fq2::set(P.valy, zero);

    let aux1 = (Yr * Zr * Zr + Yr * Zr * Zr) * Yp;
    let aux2 = -(Xr * Xr * Zr + Xr * Xr * Zr + Xr * Xr * Zr) * Xp;
    let aux3 = (Xr * Xr * Xr) + (Xr * Xr * Xr) + (Xr * Xr * Xr) - (Yr * Yr * Zr) - (Yr * Yr * Zr);

    let zero2 = Fq2::zero();

    let temp1 = Fq6::set(zero2, aux1, zero2);
    let temp2 = Fq6::set(aux3, aux2, zero2);

    Fq12::set(temp2, temp1)
}

fn ate_line_adding(RZ: G2Z, QZ: G2Z, P: G1) -> Fq12 {
    /*

    Computing L_R,Q (P)

    L: (Xq*Zr - Xr*Zq) * Yp + (YrZq - YqZr) * Xp + (XrYq - XqYr) ->

    (Xq*Zr - Xr*Zq) * const * Yp * v^2 + (YrZq - YqZr) * const * Xp * v * w +
    (XrYq - XqYr) * const * w,

    where const = (1-u)/2 = sextic_const <--- no need to multiply due to Frobenius endomorphism

    */
    let Xq = QZ.valX;
    let Yq = QZ.valY;
    let Zq = QZ.valZ;

    let Xr = RZ.valX;
    let Yr = RZ.valY;
    let Zr = RZ.valZ;

    let zero = Fq::zero();

    let Xp = Fq2::set(P.valx, zero);
    let Yp = Fq2::set(P.valy, zero);

    let aux1 = (Xq * Zr - Xr * Zq) * Yp;
    let aux2 = (Yr * Zq - Yq * Zr) * Xp;
    let aux3 = (Xr * Yq - Xq * Yr);

    let zero2 = Fq2::zero();

    let temp1 = Fq6::set(zero2, zero2, aux1);
    let temp2 = Fq6::set(aux3, aux2, zero2);

    Fq12::set(temp1, temp2)
}

//Miller's loop
fn ate(P: G1, QZ: G2Z) -> Fq12 {
    let zero = Fq::zero();
    let mut mask: u32 =
        P.inf | (Fq::is_equal(QZ.valZ.val0, zero) & Fq::is_equal(QZ.valZ.val1, zero));
    mask = mask.wrapping_neg(); //0xF...F if either P or Q is inf; 0, otherwise

    let mut RZ = QZ;

    let mut f = Fq12::one().mont_rep();
    let mut L = Fq12::zero();

    for bit in bATE_STR.chars().skip(1) {
        L = ate_line_doubling(RZ, P);
        RZ = RZ + RZ;
        f = f * f;
        f = f * L;
        if bit == '1' {
            L = ate_line_adding(RZ, QZ, P);
            RZ = RZ + QZ;
            f = f * L;
        }
    }

    f.val0.val0.val0.limbs[0] = (f.val0.val0.val0.limbs[0] & !mask) | (1u32 & mask);
    f
}

pub fn ate_pairing(P: G1, QZ: G2Z) -> Fq12 {
    frobenius::ate_exp(ate(P, QZ))
}
