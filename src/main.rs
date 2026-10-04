use nalgebra::DVector;
use num_complex::Complex;

type StateVector = DVector<Complex<f64>>;

fn main() {
    let n = 8;
    let mut state: StateVector = DVector::from_element(n, Complex::new(0.0, 0.0));
    state[0] = Complex::new(1.0, 0.0);

    println!("state = {:?}", state);
}
