package fields.polynomials;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

class PolynomialLinearEquationsTest {

	@Test
	void testBasics() throws IOException {
		PrimeField f257 = PrimeField.getPrimeField(257);
		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(f257, 2, Monomial.GREVLEX);
		List<List<Polynomial<PFE>>> matrixList = new ArrayList<>();
		List<Polynomial<PFE>> row1 = new ArrayList<>();
		row1.add(polynomials.parse("X*Y"));
		row1.add(polynomials.parse("X^2"));
		row1.add(polynomials.parse("Y^2"));
		matrixList.add(row1);
		List<Polynomial<PFE>> row2 = new ArrayList<>();
		row2.add(polynomials.parse("Y^2 + -1*X^3 + X"));
		row2.add(polynomials.parse("X*Y"));
		row2.add(polynomials.parse("0"));
		matrixList.add(row2);
		Matrix<Polynomial<PFE>> matrix = new Matrix<>(matrixList);
		MatrixModule<Polynomial<PFE>> mm = matrix.getModule(polynomials);
		List<Vector<Polynomial<PFE>>> syzygies = polynomials.syzygyProblem(mm, matrix);
		for (Vector<Polynomial<PFE>> syzygy : syzygies) {
			assertEquals(mm.codomain().zero(), mm.multiply(matrix, syzygy));
		}
	}

	@Test
	void testZeroRow() throws IOException {
		PrimeField f257 = PrimeField.getPrimeField(257);
		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(f257, 2, Monomial.GREVLEX);
		List<List<Polynomial<PFE>>> matrixList = new ArrayList<>();
		List<Polynomial<PFE>> row1 = new ArrayList<>();
		row1.add(polynomials.parse("X*Y"));
		row1.add(polynomials.parse("X^2"));
		row1.add(polynomials.parse("Y^2"));
		matrixList.add(row1);
		List<Polynomial<PFE>> row2 = new ArrayList<>();
		row2.add(polynomials.parse("0"));
		row2.add(polynomials.parse("0"));
		row2.add(polynomials.parse("0"));
		matrixList.add(row2);
		Matrix<Polynomial<PFE>> matrix = new Matrix<>(matrixList);
		MatrixModule<Polynomial<PFE>> mm = matrix.getModule(polynomials);
		List<Vector<Polynomial<PFE>>> syzygies = polynomials.syzygyProblem(mm, matrix);
		System.out.println(syzygies);
		for (Vector<Polynomial<PFE>> syzygy : syzygies) {
			System.out.println(mm.multiply(matrix, syzygy));
		}
	}

}
