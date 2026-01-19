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

use std::sync::LazyLock;

const bX_STR: &str = "1101001000000001000000000000000000000000000000010000000000000000"; // |x|
const bXP1_STR: &str = "1101001000000001000000000000000000000000000000010000000000000001"; // |x| + 1
const bXP1_DIV3_STR: &str = "100011000000000010101010101010101010101010101011010101010101011"; // (1 + |x|)/3
const bX2_STR: &str = "10101100010001011010010000000001000000000000000110100100000000100000000000000000000000000000000100000000000000000000000000000000"; //x^2

/*

    If f = h + gw, then f = g0 + h0 * w + g1 * w^2 + h1 * w^3 + g2 * w^4 + h2 * w^5.

    f^q = ~g0 + ~h0 * g11 * w + ~g1 * g12 *  w^2 + ~h1 * g13 * w^3 + ~g2 * g14 * w^4 + ~h2 * g15 * w^5, where ~a is conjugation of a.

    f^(q^2) = f = g0 + h0 * g21 * w + g1 * g22 * w^2 + h1 * g23 * w^3 + g2 * g24 * w^4 + h2 * g25 * w^5.

    Coefficients g11, ..., g15, g21, ..., g25 are given below.

    Notice that g2j = g1j * ~g1j.

    Fq2_set_hex_str(&g11, "1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8",
    "fc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3");

    Fq2_set_hex_str(&g12, "0", "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac");

    Fq2_set_hex_str(&g13, "6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09",
    "6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09");

    Fq2_set_hex_str(&g14, "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad", "0");

    Fq2_set_hex_str(&g15,  "5b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116",
    "144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995");

    Fq2_set_hex_str(&g21, "5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff", "0");
    Fq2_set_hex_str(&g22, "5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe", "0");
    Fq2_set_hex_str(&g23, "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa", "0");
    Fq2_set_hex_str(&g24, "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac", "0");
    Fq2_set_hex_str(&g25, "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad", "0");

*/

const g11: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0x92235fb8, 0x8d0775ed, 0x63e7813d, 0xf67ea53d, 0x84bab9c4, 0x7b2443d7, 0x3cbd5f4f,
                0x0fd603fd, 0x202c0d1f, 0xc231beb4, 0x02bb0667, 0x1904d3bf,
            ],
        },
        val1: Fq {
            limbs: [
                0x6ddc4af3, 0x2cf78a12, 0x4d6c7ec2, 0x282d5ac1, 0x71f63c5f, 0xec0c8ec9, 0xb6c7b36f,
                0x54a14787, 0x231f9fb8, 0x88e9e902, 0x36c4e032, 0x00fc3e2b,
            ],
        },
    }
};

const g12: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [0u32; N_LIMBS],
        },

        val1: Fq {
            limbs: [
                0x0000aaac, 0x8bfd0000, 0x4f49fffd, 0x409427eb, 0x0fb85f9b, 0x897d2965, 0x89759ad4,
                0xaa0d857d, 0x63d4de85, 0xec024086, 0x397fe699, 0x1a0111ea,
            ],
        },
    }
};

const g13: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0xede3cc09, 0xc81084fb, 0x72ec05f4, 0xee67992f, 0x009241c5, 0x77f76e17, 0xc2d3435e,
                0x48395dab, 0x6bd17ffe, 0x6831e36d, 0x37ff400b, 0x06af0e04,
            ],
        },
        val1: Fq {
            limbs: [
                0xede3cc09, 0xc81084fb, 0x72ec05f4, 0xee67992f, 0x009241c5, 0x77f76e17, 0xc2d3435e,
                0x48395dab, 0x6bd17ffe, 0x6831e36d, 0x37ff400b, 0x06af0e04,
            ],
        },
    }
};

const g14: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0x0000aaad, 0x8bfd0000, 0x4f49fffd, 0x409427eb, 0x0fb85f9b, 0x897d2965, 0x89759ad4,
                0xaa0d857d, 0x63d4de85, 0xec024086, 0x397fe699, 0x1a0111ea,
            ],
        },

        val1: Fq {
            limbs: [0u32; N_LIMBS],
        },
    }
};

const g15: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0x80078116, 0x9b18fae9, 0x257f8732, 0xc63a3e6e, 0x8e9c0566, 0x8beadf4d, 0x0c0b8fee,
                0xf3981624, 0x48b1e045, 0xdf47fa6b, 0x013a5fd8, 0x05b2cfd9,
            ],
        },

        val1: Fq {
            limbs: [
                0x7ff82995, 0x1ee60516, 0x8bd478cd, 0x5871c190, 0x6814f0bd, 0xdb45f353, 0xe77982d0,
                0x70df3560, 0xfa99cc91, 0x6bd3ad4a, 0x384586c1, 0x144e4211,
            ],
        },
    }
};

const g21: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0xfffeffff, 0x2e01ffff, 0x620a0002, 0xde17d813, 0xe6f89688, 0xddb3a93b, 0x6a0f77ea,
                0xba69c607, 0xdf76ce51, 0x5f19672f, 0x00000000, 0x00000000,
            ],
        },

        val1: Fq {
            limbs: [0u32; N_LIMBS],
        },
    }
};

const g22: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0xfffefffe, 0x2e01ffff, 0x620a0002, 0xde17d813, 0xe6f89688, 0xddb3a93b, 0x6a0f77ea,
                0xba69c607, 0xdf76ce51, 0x5f19672f, 0x00000000, 0x00000000,
            ],
        },

        val1: Fq {
            limbs: [0u32; N_LIMBS],
        },
    }
};

const g23: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0xffffaaaa, 0xb9feffff, 0xb153ffff, 0x1eabfffe, 0xf6b0f624, 0x6730d2a0, 0xf38512bf,
                0x64774b84, 0x434bacd7, 0x4b1ba7b6, 0x397fe69a, 0x1a0111ea,
            ],
        },

        val1: Fq {
            limbs: [0u32; N_LIMBS],
        },
    }
};

const g24: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0x0000aaac, 0x8bfd0000, 0x4f49fffd, 0x409427eb, 0x0fb85f9b, 0x897d2965, 0x89759ad4,
                0xaa0d857d, 0x63d4de85, 0xec024086, 0x397fe699, 0x1a0111ea,
            ],
        },

        val1: Fq {
            limbs: [0u32; N_LIMBS],
        },
    }
};

const g25: Fq2 = {
    Fq2 {
        val0: Fq {
            limbs: [
                0x0000aaad, 0x8bfd0000, 0x4f49fffd, 0x409427eb, 0x0fb85f9b, 0x897d2965, 0x89759ad4,
                0xaa0d857d, 0x63d4de85, 0xec024086, 0x397fe699, 0x1a0111ea,
            ],
        },

        val1: Fq {
            limbs: [0u32; N_LIMBS],
        },
    }
};

fn ate_exp_q(a: Fq12) -> Fq12 {
    /*

    f^q = ~g0 + ~h0 * g11 * w + ~g1 * g12 *  w^2 + ~h1 * g13 * w^3 + ~g2 * g14 * w^4 + ~h2 * g15 * w^5

    */

    let g0 = a.val0.val0.conj();
    let g1 = a.val0.val1.conj() * g12.mont_rep();
    let g2 = a.val0.val2.conj() * g14.mont_rep();

    let h0 = a.val1.val0.conj() * g11.mont_rep();
    let h1 = a.val1.val1.conj() * g13.mont_rep();
    let h2 = a.val1.val2.conj() * g15.mont_rep();

    Fq12::set(Fq6::set(g0, g1, g2), Fq6::set(h0, h1, h2))
}

fn ate_exp_q2(a: Fq12) -> Fq12 {
    /*

    f^(q^2) = g0 + h0 * g21 * w + g1 * g22 * w^2 + h1 * g23 * w^3 + g2 * g24 * w^4 + h2 * g25 * w^5.

    */

    let g0 = a.val0.val0;
    let g1 = a.val0.val1 * g22.mont_rep();
    let g2 = a.val0.val2 * g24.mont_rep();

    let h0 = a.val1.val0 * g21.mont_rep();
    let h1 = a.val1.val1 * g23.mont_rep();
    let h2 = a.val1.val2 * g25.mont_rep();

    Fq12::set(Fq6::set(g0, g1, g2), Fq6::set(h0, h1, h2))
}

fn ate_exp_q6(a: Fq12) -> Fq12 {
    /*

    f = g + hw

    f^(q^6) = g - hw

    */
    Fq12::set(a.val0, -a.val1)
}

fn ate_exp_x(a: Fq12) -> Fq12 {
    let mut aux = a;

    for bit in bX_STR.chars().skip(1) {
        aux = aux * aux;
        if bit == '1' {
            aux = aux * a;
        }
    }

    aux
}

fn ate_exp_xp1(a: Fq12) -> Fq12 {
    let mut aux = a;

    for bit in bXP1_STR.chars().skip(1) {
        aux = aux * aux;
        if bit == '1' {
            aux = aux * a;
        }
    }

    aux
}

fn ate_exp_xp1div3(a: Fq12) -> Fq12 {
    let mut aux = a;

    for bit in bXP1_DIV3_STR.chars().skip(1) {
        aux = aux * aux;
        if bit == '1' {
            aux = aux * a;
        }
    }

    aux
}

//rop = f^((q^4 - q^2 + 1)/r)
fn ate_exp_hard(f: Fq12) -> Fq12 {
    /*

    Note the following identity.

    (q^4 - q^2 + 1)/r = (x+1)^2 /3 * (q-x)(x^2 + q^2 - 1) + 1, where x = 0xd201000000010000 (x is taken positive here).

    Also, f^(-1) = ~f, because of the consequent raising to the power (q^6-1).

    */

    // I. powering to (x+1)^2/3
    let a = ate_exp_xp1(ate_exp_xp1div3(f));
    // II. powering to q - x
    let b = ate_exp_q(a) * ate_exp_x(a.conj());
    // III. powering to x^2 + q^2 - 1
    let c = ate_exp_x(ate_exp_x(b)) * ate_exp_q2(b) * b.conj();
    // IV. multiplying by f and returning the result
    c * f
}

pub fn ate_exp(f: Fq12) -> Fq12 {
    // Computing f^(q^2 + 1)
    let a = f * ate_exp_q2(f);
    // Computing (a^(q^2 + 1)) ^ (q^6 - 1)
    let b = ate_exp_q6(a) * a.mont_mul_inv();
    // Powering to the (q^4 - q^2 + 1)/r
    let c = ate_exp_hard(b);
    c
}
