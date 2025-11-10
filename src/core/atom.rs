pub type AtomId = usize;

macro_rules! define_elements {
    ($($name:ident = $value:literal),* $(,)?) => {
        #[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
        #[repr(u8)]
        pub enum Element {
            $($name = $value),*
        }

        impl Element {
            pub fn from_atomic_number(atomic_number: u8) -> Option<Self> {
                match atomic_number {
                    $($value => Some(Element::$name),)*
                    _ => None,
                }
            }

            pub fn atomic_number(self) -> u8 {
                self as u8
            }

            pub fn valence_electrons(self) -> Option<u8> {
                use Element::*;
                let electrons = match self {
                    H => 1, He => 2,
                    Li | Na | K | Rb | Cs | Fr => 1,
                    Be | Mg | Ca | Sr | Ba | Ra => 2,
                    B | Al | Ga | In | Tl => 3,
                    C | Si | Ge | Sn | Pb => 4,
                    N | P | As | Sb | Bi => 5,
                    O | S | Se | Te | Po => 6,
                    F | Cl | Br | I | At => 7,
                    Ne | Ar | Kr | Xe | Rn | Og => 8,
                    _ => return None,
                };
                Some(electrons)
            }

            pub fn is_common_conjugation_element(self) -> bool {
                matches!(
                    self,
                    Element::B | Element::C | Element::N | Element::O | Element::P |
                    Element::S | Element::Si | Element::As | Element::Se | Element::Ge |
                    Element::F | Element::Cl | Element::Br | Element::I
                )
            }
        }

        impl std::str::FromStr for Element {
            type Err = String;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                let s = s.trim();
                // Try parse as atomic number first
                if let Ok(n) = s.parse::<u8>() {
                    if let Some(e) = Element::from_atomic_number(n) {
                        return Ok(e);
                    }
                }
                // Try match symbol (case-insensitive)
                $(
                    if s.eq_ignore_ascii_case(stringify!($name)) {
                        return Ok(Element::$name);
                    }
                )*
                Err(format!("invalid element: {}", s))
            }
        }
    };
}

define_elements!(
    H = 1,
    He = 2,
    Li = 3,
    Be = 4,
    B = 5,
    C = 6,
    N = 7,
    O = 8,
    F = 9,
    Ne = 10,
    Na = 11,
    Mg = 12,
    Al = 13,
    Si = 14,
    P = 15,
    S = 16,
    Cl = 17,
    Ar = 18,
    K = 19,
    Ca = 20,
    Sc = 21,
    Ti = 22,
    V = 23,
    Cr = 24,
    Mn = 25,
    Fe = 26,
    Co = 27,
    Ni = 28,
    Cu = 29,
    Zn = 30,
    Ga = 31,
    Ge = 32,
    As = 33,
    Se = 34,
    Br = 35,
    Kr = 36,
    Rb = 37,
    Sr = 38,
    Y = 39,
    Zr = 40,
    Nb = 41,
    Mo = 42,
    Tc = 43,
    Ru = 44,
    Rh = 45,
    Pd = 46,
    Ag = 47,
    Cd = 48,
    In = 49,
    Sn = 50,
    Sb = 51,
    Te = 52,
    I = 53,
    Xe = 54,
    Cs = 55,
    Ba = 56,
    La = 57,
    Ce = 58,
    Pr = 59,
    Nd = 60,
    Pm = 61,
    Sm = 62,
    Eu = 63,
    Gd = 64,
    Tb = 65,
    Dy = 66,
    Ho = 67,
    Er = 68,
    Tm = 69,
    Yb = 70,
    Lu = 71,
    Hf = 72,
    Ta = 73,
    W = 74,
    Re = 75,
    Os = 76,
    Ir = 77,
    Pt = 78,
    Au = 79,
    Hg = 80,
    Tl = 81,
    Pb = 82,
    Bi = 83,
    Po = 84,
    At = 85,
    Rn = 86,
    Fr = 87,
    Ra = 88,
    Ac = 89,
    Th = 90,
    Pa = 91,
    U = 92,
    Np = 93,
    Pu = 94,
    Am = 95,
    Cm = 96,
    Bk = 97,
    Cf = 98,
    Es = 99,
    Fm = 100,
    Md = 101,
    No = 102,
    Lr = 103,
    Rf = 104,
    Db = 105,
    Sg = 106,
    Bh = 107,
    Hs = 108,
    Mt = 109,
    Ds = 110,
    Rg = 111,
    Cn = 112,
    Nh = 113,
    Fl = 114,
    Mc = 115,
    Lv = 116,
    Ts = 117,
    Og = 118,
);
