package fields.vectors;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.FiniteComplexVectorSpace;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.interfaces.InnerProductSpace.SesquilinearFormDiagonalizationResult;
import fields.interfaces.SesquilinearForm;

class SesquilinearFormTest {

	@Test
	void testComplex() {
		Complex c = Complex.c(128);
		int n = 5;
		List<List<ComplexNumber>> coeffs = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			List<ComplexNumber> row = new ArrayList<>();
			for (int j = 0; j < i; j++) {
				row.add(c.conjugate(coeffs.get(j).get(i)));
			}
			row.add(c.getEmbedding(c.getReals().getRandomElement()));
			for (int j = i + 1; j < n; j++) {
				row.add(c.getRandomElement());
			}
			coeffs.add(row);
		}
		Matrix<ComplexNumber> matrix = new Matrix<>(coeffs);
		FiniteComplexVectorSpace space = new FiniteComplexVectorSpace(c, n);
		SesquilinearForm<ComplexNumber, Vector<ComplexNumber>> form = space.getSesquilinearForm(matrix);
		SesquilinearFormDiagonalizationResult<ComplexNumber, Vector<ComplexNumber>> diag = space
				.diagonalizeSesquilinearForm(matrix);
		for (int i = 0; i < diag.getPositive(); i++) {
			System.out.println(form.evaluate(diag.getBasis().get(i), diag.getBasis().get(i)) + " =  1");
		}
		for (int i = diag.getPositive(); i < diag.getPositive() + diag.getNegative(); i++) {
			System.out.println(form.evaluate(diag.getBasis().get(i), diag.getBasis().get(i)) + " = -1");
		}
		for (int i = diag.getPositive() + diag.getNegative(); i < n; i++) {
			System.out.println(form.evaluate(diag.getBasis().get(i), diag.getBasis().get(i)) + " =  0");
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					continue;
				}
				System.out.println(form.evaluate(diag.getBasis().get(i), diag.getBasis().get(j)) + " = 0");
			}
		}
	}

	@Test
	void testReal() {
		Reals r = Reals.r(128);
		int n = 5;
		List<List<Real>> coeffs = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			List<Real> row = new ArrayList<>();
			for (int j = 0; j < i; j++) {
				row.add(coeffs.get(j).get(i));
			}
			for (int j = i; j < n; j++) {
				row.add(r.getRandomElement());
			}
			coeffs.add(row);
		}
		Matrix<Real> matrix = new Matrix<>(coeffs);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, n);
		SesquilinearFormDiagonalizationResult<Real, Vector<Real>> diag = space.diagonalizeSesquilinearForm(matrix);
	}

}
