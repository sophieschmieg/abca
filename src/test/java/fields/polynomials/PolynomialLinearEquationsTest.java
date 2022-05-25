package fields.polynomials;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

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
		Vector<Polynomial<PFE>> rhs = new Vector<>(polynomials.parse("1"), polynomials.parse("0"));
		assertFalse(polynomials.isSubModuleMember(mm, matrix, rhs));
		rhs = mm.multiply(matrix,  new Vector<>(polynomials.parse("Y"), polynomials.parse("X^2 + -1"), polynomials.parse("-1*X")));
		assertTrue(polynomials.isSubModuleMember(mm, matrix, rhs));
		assertEquals(rhs, mm.multiply(matrix, polynomials.asSubModuleMember(mm, matrix, rhs)));
	}

	@Test
	void testZeroRow2() throws IOException {
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
		for (Vector<Polynomial<PFE>> syzygy : syzygies) {
			assertEquals(mm.codomain().zero(), mm.multiply(matrix, syzygy));
		}
		Vector<Polynomial<PFE>> rhs = new Vector<>(polynomials.parse("1"), polynomials.parse("0"));
		assertFalse(polynomials.isSubModuleMember(mm, matrix, rhs));
		rhs = new Vector<>(polynomials.parse("0"), polynomials.parse("Y^2 + 3*X^2*Y"));
		assertFalse(polynomials.isSubModuleMember(mm, matrix, rhs));
		rhs = new Vector<>( polynomials.parse("Y^2 + 3*X^2*Y"),polynomials.parse("0"));
		assertTrue(polynomials.isSubModuleMember(mm, matrix, rhs));
			assertEquals(rhs, mm.multiply(matrix, polynomials.asSubModuleMember(mm, matrix, rhs)));
	
	}

	@Test
	void testZeroRow1() throws IOException {
		PrimeField f257 = PrimeField.getPrimeField(257);
		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(f257, 2, Monomial.GREVLEX);
		List<List<Polynomial<PFE>>> matrixList = new ArrayList<>();
		List<Polynomial<PFE>> row1 = new ArrayList<>();
		row1.add(polynomials.parse("0"));
		row1.add(polynomials.parse("0"));
		row1.add(polynomials.parse("0"));
		matrixList.add(row1);
		List<Polynomial<PFE>> row2 = new ArrayList<>();
		row2.add(polynomials.parse("X*Y"));
		row2.add(polynomials.parse("X^2"));
		row2.add(polynomials.parse("Y^2"));
		matrixList.add(row2);
		Matrix<Polynomial<PFE>> matrix = new Matrix<>(matrixList);
		MatrixModule<Polynomial<PFE>> mm = matrix.getModule(polynomials);
		List<Vector<Polynomial<PFE>>> syzygies = polynomials.syzygyProblem(mm, matrix);
		for (Vector<Polynomial<PFE>> syzygy : syzygies) {
			assertEquals(mm.codomain().zero(), mm.multiply(matrix, syzygy));
		}
		Vector<Polynomial<PFE>> rhs = new Vector<>(polynomials.parse("1"), polynomials.parse("0"));
		assertFalse(polynomials.isSubModuleMember(mm, matrix, rhs));
		rhs = new Vector<>(polynomials.parse("0"), polynomials.parse("Y^2 + 3*X^2*Y"));
		assertTrue(polynomials.isSubModuleMember(mm, matrix, rhs));
		assertEquals(rhs, mm.multiply(matrix, polynomials.asSubModuleMember(mm, matrix, rhs)));
	}
}
