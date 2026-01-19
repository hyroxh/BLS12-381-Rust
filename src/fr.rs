use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};
use std::ptr;

pub const LEN_MAX_R: usize = 64;
pub const FR_LIMBS: usize = 8;

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Fr {
    pub limbs: [u32; FR_LIMBS],
}

pub const R_MODULUS: Fr = Fr {
    limbs: [
        0x00000001, 0xffffffff, 0xfffe5bfe, 0x53bda402, 0x09a1d805, 0x3339d808, 0x299d7d48,
        0x73eda753,
    ],
};

impl Fr {
	
    // Parse a hex string into Fr
    pub fn from_hex_str(hex_str: &str) -> Result<Self, &'static str> {
        let mut limbs = [0u32; FR_LIMBS];
        let clean_str = hex_str.trim_start_matches("0x");

        // We process the string from right to left (LSB to MSB)
        let chars: Vec<char> = clean_str.chars().collect();
        let len = chars.len();

        for i in 0..FR_LIMBS {
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

        Ok(Fr { limbs })
    }
}

// --- Display / Printing ---

impl fmt::Display for Fr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "0x")?;
        // Print from Most Significant Limb to Least
        for i in (0..FR_LIMBS).rev() {
            write!(f, "{:08x}", self.limbs[i])?;
        }
        Ok(())
    }
}
