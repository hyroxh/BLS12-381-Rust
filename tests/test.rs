use bls12_381_rust::fq::Fq;
use bls12_381_rust::fq2::Fq2;
use bls12_381_rust::fq6::Fq6;
use bls12_381_rust::fq12::Fq12;
use bls12_381_rust::fr::Fr;
use bls12_381_rust::fr::R_MODULUS;
use bls12_381_rust::g1::{G1, G1Z, genG1, orderG1};
use bls12_381_rust::g2::{G2, G2Z, genG2, orderG2};
use bls12_381_rust::pairing::ate_pairing;

fn ate_pairing_mul(a: &Fr, b: &Fr) -> Fq12 {
    let P = (*a * genG1.to_project().mont_rep()).to_affine();
    let QZ = *b * genG2.to_project().mont_rep();
    ate_pairing(P, QZ).mont_rep_inv()
}

fn ate_test(a: &Fr, b: &Fr, print: u8) -> u8 {
    let res1 = ate_pairing_mul(a, b);
    let res2 = ate_pairing_mul(b, a);

    if Fq12::is_equal(res1, res2) == 1u32 {
        if print == 1 {
            println!("{}", res1);
        }
        1u8
    } else {
        0u8
    }
}

fn main() {
    //----Test1
    let a = Fr::from_hex_str("1aaaaaaaaaaffffffffbf").unwrap();
    let b = Fr::from_hex_str("12345ccccddf").unwrap();

    if ate_test(&a, &b, 1) == 1 {
        println!("Test 1:\tSuccess!");
        println!("----------");
    }

    //----Test2
    let a = Fr::from_hex_str("1").unwrap();
    let b = R_MODULUS;

    if ate_test(&a, &b, 1) == 1 {
        println!("Test 2:\tSuccess!");
        println!("----------");
    }

    //----Test3
    let a = Fr::from_hex_str("aaaaaaaaaa33333333333333333333").unwrap();
    let b = Fr::from_hex_str("1").unwrap();

    if ate_test(&a, &b, 0) == 1 {
        println!("Test 3:\tSuccess!");
        println!("----------");
    }

    //----Test3
    let a = Fr::from_hex_str("1").unwrap();
    let b = Fr::from_hex_str("ffffffffffffffffffff7777").unwrap();

    if ate_test(&a, &b, 0) == 1 {
        println!("Test 4:\tSuccess!");
        println!("----------");
    }
}
