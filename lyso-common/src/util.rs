use std::io::BufRead;

pub static DNA: [char; 5] = ['A', 'T', 'G', 'C', 'N'];

pub fn is_dna(c: char) -> bool {
    matches!(c, 'A' | 'T' | 'G' | 'C' | 'N')
}

// TODO accept an arbitrary number of validator functions
pub trait Validate {
    fn valid(&self) -> Result<bool, &'static str> {
        let sv = self.seq_valid();
        let svb = match sv {
            Ok(val) => val,
            Err(e) => return Err(e),
        };

        let qv = self.qual_valid();
        let qvb = match qv {
            Ok(val) => val,
            Err(e) => return Err(e),
        };

        Ok(svb && qvb)
    }

    fn seq_valid(&self) -> Result<bool, &'static str>;

    fn qual_valid(&self) -> Result<bool, &'static str>;
}


// --- BEGIN TESTS --- //

#[cfg(test)]
mod tests {
    use crate::util;
    #[test]
    fn valid_dna() {
        // We only test uppercase because all strings
        // will be capitalized a-priori
        let good_dna = "ATGCNCTGA";
        let bad_dna = "XYZLKJM1234";
        for c in good_dna.chars() {
            assert!(util::is_dna(c));
        }

        for c in bad_dna.chars() {
            assert_ne!(util::is_dna(c), true);
        }
    }
}

// --- END TESTS --- //
