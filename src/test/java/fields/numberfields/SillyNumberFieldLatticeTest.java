package fields.numberfields;

import static org.junit.jupiter.api.Assertions.fail;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractLattice;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Lattice;
import fields.interfaces.MathMap;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField.NFE;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.SubRealInnerProductSpace;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.MiscAlgorithms;

class SillyNumberFieldLatticeTest {

	private MathMap<NFE, NFE> automorphism(int i, NumberField field, List<NFE> conjugates) {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomialRing = q.getUnivariatePolynomialRing();
		return new FunctionMathMap<>((NFE t) -> {
			UnivariatePolynomial<Fraction> asPolynomial = t.asPolynomial();
			asPolynomial = polynomialRing.substitute(asPolynomial,
					Collections.singletonList(conjugates.get(i).asPolynomial()));
			return field.fromPolynomial(asPolynomial);
		});
	}

	@Test
	void test() {
		int n = 4;
		int k = 3;
		int q = 41;
		Integers z = Integers.z();
		Rationals rat = Rationals.q();
		// 2cos(2phi) = 2cos(phi)^2 - 2sin(phi)^2
		// = 4cos(phi)^2 - 2 = (2cos(phi))^2 - 2
		// 2cos(pi/2) = 0
		UnivariatePolynomialRing<IntE> polynomialRing = z.getUnivariatePolynomialRing();
		UnivariatePolynomial<IntE> twoCos = polynomialRing.getVar();
		for (int i = 0; i < n - 2; i++) {
			twoCos = polynomialRing.subtract(polynomialRing.power(twoCos, 2), polynomialRing.getInteger(2));
		}
		NumberField field = NumberField.getNumberFieldFromIntegerPolynomial(twoCos);
		// x^(2^n) + 1 = 0
		// 2cos(phi) = x + x^-1
		// sigma(x) = x^k for k 2i+1, i = 0..2^n
		// sigma(x+x^-1) = x^k + x^-k
		// (x+x^-1)^k = sum binom(k, j)*x^j*x^-(k-j) = sum binom(k, j) x^(2j-k), j=0..k
		// sum binom(k, j)sigma_j(2cos(phi))
		List<NFE> conjugates = new ArrayList<>();
		conjugates.add(field.alpha());
		for (int i = 1; i < field.degree(); i++) {
			NFE conjugate = field.power(field.alpha(), 2 * i + 1);
			for (int j = 0; j < i; j++) {
				conjugate = field.subtract(conjugate,
						field.multiply(MiscAlgorithms.binomial(2 * i + 1, i - j), conjugates.get(j)));
			}
			conjugates.add(conjugate);
		}

		System.out.println(field);
		NFE unit = field.add(field.alpha(), field.one());
		List<NFE> units = new ArrayList<>();
		List<Vector<Real>> logs = new ArrayList<>();
		for (int i = 0; i < field.degree() - 1; i++) {
			units.add(automorphism(i, field, conjugates).evaluate(unit));
			logs.add(field.logRepresentation(units.get(i)));
		}
		Lattice<NFE, Real, Vector<Real>> unitLattice = new AbstractLattice<NFE, Real, Vector<Real>>(
				field.logRepresentationSpace(), units, true) {

			@Override
			public Vector<Real> embedding(NFE t) {
				return field.logRepresentation(t);
			}

			@Override
			public NFE zero() {
				return field.one();
			}

			@Override
			public NFE add(NFE s1, NFE s2) {
				return field.multiply(s1, s2);
			}

			@Override
			public NFE negative(NFE s) {
				return field.inverse(s);
			}

			@Override
			public NFE scalarMultiply(IntE t, NFE s) {
				return field.power(s, t);
			}

			@Override
			public Exactness exactness() {
				return field.exactness();
			}
		};
		NumberFieldIntegers order = field.maximalOrder();
		for (NFE t : unitLattice.getModuleGenerators()) {
			System.out.println(t);
			System.out.println(order.isElement(t));
			System.out.println(order.isUnit(t));
		}
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				for (int m = 0; m < 10; m++) {
					for (int l = 0; l < 10; l++) {
						if (i == 0 && j == 0 && m == 0 && l == 0) {
							continue;
						}
						Vector<IntE> testVector = new Vector<>(z.getInteger(i), z.getInteger(j), z.getInteger(m),
								z.getInteger(l));
						NFE t = field.inverse(order.fromVector(testVector));
						Vector<Real> log = field.logRepresentation(t);
						NFE u = field.logRepresentationSpace().closestLatticePoint(log, unitLattice);
						NFE s = field.divide(t, u);
						NFE rounded = field.minkowskiEmbeddingSpace().closestLatticePoint(field.minkowskiEmbedding(s),
								order);
						Fraction result = field.norm(field.subtract(t, rounded));
						if (result.compareTo(Rationals.q().one()) > 0) {
							System.err.println("Rounded " + t + " (" + i + ", " + j + ", " + m + ", " + l + ", "
									+ ") to " + rounded + " but norm is " + result);
						}
					}
				}
			}
		}
		fail();
		PrimeField fp = PrimeField.getPrimeField(q);
		UnivariatePolynomialRing<PFE> reducedPolynomialRing = fp.getUnivariatePolynomialRing();
		List<List<NFE>> matrixList = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			List<NFE> row = new ArrayList<>();
			for (int j = 0; j < k; j++) {
				UnivariatePolynomial<IntE> value = z
						.liftUnivariatePolynomial(reducedPolynomialRing.getRandomElement(n));
				row.add(field
						.fromPolynomial(rat.getUnivariatePolynomialRing().getEmbedding(value, rat.getEmbeddingMap())));
			}
			for (int j = 0; j < k; j++) {
				row.add(i == j ? field.one() : field.zero());
			}
			for (int j = 0; j < k; j++) {
				row.add(i == j ? field.getInteger(q) : field.zero());
			}
			matrixList.add(row);
		}
		Matrix<NFE> matrix = new Matrix<>(matrixList);
		List<Vector<NFE>> kernel = order.syzygyProblem(matrix);
		List<Vector<NFE>> projectedKernel = new ArrayList<>();
		for (Vector<NFE> v : kernel) {
			projectedKernel.add(new Vector<>(v.asList().subList(0, 2 * k)));
		}
		System.out.println(projectedKernel.size());
		System.out.println(projectedKernel);
		List<Vector<NFE>> reduced = computeLatticeReduction(projectedKernel, field, new FreeModule<>(field, 2 * k),
				unitLattice);
		System.out.println(reduced);
	}

	private List<Vector<NFE>> computeLatticeReduction(List<Vector<NFE>> sublatticeBasis, NumberField field,
			FreeModule<NFE> module, Lattice<NFE, Real, Vector<Real>> unitLattice) {
		Integers z = Integers.z();
		SubRealInnerProductSpace<Real, Vector<Real>> projection = new SubRealInnerProductSpace<>(
				field.logRepresentationSpace(), unitLattice.generatorsAsMatrix().asColumnList());
		List<Vector<NFE>> basis = new ArrayList<>();
		basis.addAll(sublatticeBasis);
		List<Vector<NFE>> orthogonal = new ArrayList<>();
		List<NFE> values = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			values.add(field.zero());
			orthogonal.add(module.zero());
		}
		List<List<NFE>> coefficients = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			coefficients.add(new ArrayList<>());
			for (int j = 0; j < i; j++) {
				coefficients.get(i).add(field.zero());
			}
		}
		int k = 1;
		int kmax = 0;
		orthogonal.set(0, basis.get(0));
		values.set(0, module.innerProduct(orthogonal.get(0), orthogonal.get(0)));
		while (k < basis.size()) {
			if (values.get(k - 1).equals(field.zero())) {
				basis.remove(k - 1);
				coefficients.remove(coefficients.size() - 1);
				orthogonal.remove(k - 1);
				values.remove(k - 1);
				k = Math.max(1, k - 1);
				kmax = Math.min(kmax, k - 1);
				continue;
			}
			if (k > kmax) {
				kmax = k;
				orthogonal.set(k, basis.get(k));
				Vector<NFE> embedded = orthogonal.get(k);
				for (int j = 0; j < k; j++) {
					coefficients.get(k).set(j,
							field.divide(module.innerProduct(embedded, orthogonal.get(j)), values.get(j)));
					orthogonal.set(k, module.subtract(orthogonal.get(k),
							module.scalarMultiply(coefficients.get(k).get(j), orthogonal.get(j))));
				}
				values.set(k, module.innerProduct(orthogonal.get(k), orthogonal.get(k)));
			}
			reduceLLLBasis(k, k - 1, coefficients, basis, field, module);
			if (values.get(k).equals(field.zero())) {
				basis.remove(k);
				coefficients.remove(coefficients.size() - 1);
				orthogonal.remove(k);
				values.remove(k);
				k = Math.max(1, k - 1);
				kmax = Math.min(kmax, k - 1);
				continue;
			}
			FiniteRealVectorSpace space = field.minkowskiEmbeddingSpace();
			Vector<Real> kValue = field.logRepresentation(values.get(k));
			Vector<Real> adjustedValueFactor1 = space.subtract(field.minkowskiEmbedding(field.one()),
					field.minkowskiEmbedding(
							field.multiply(coefficients.get(k).get(k - 1), coefficients.get(k).get(k - 1))));
			Vector<Real> adjustedValueFactor2 = field.minkowskiEmbedding(values.get(k - 1));
			List<Real> adjustedValueList = new ArrayList<>();
			Reals r = space.getValueField();
			for (int i = 0; i < field.degree(); i++) {
				adjustedValueList.add(
						r.log(r.abs(r.multiply(adjustedValueFactor1.get(i + 1), adjustedValueFactor2.get(i + 1)))));
			}
			Vector<Real> adjustedValue = new Vector<>(adjustedValueList);
			boolean comparable = true;
			boolean inverted = false;
			boolean invertedSet = false;
			Real kSum = r.zero();
			Real adjustedSum = r.zero();
			for (int i = 0; i < field.degree(); i++) {
				adjustedSum = r.add(adjustedSum, adjustedValue.get(i + 1));
				kSum = r.add(kSum, kValue.get(i + 1));
				int cmp = kValue.get(i + 1).compareTo(adjustedValue.get(i + 1));
				if (cmp == 0) {
					continue;
				}
				if (!invertedSet) {
					invertedSet = true;
					inverted = cmp < 0;
					continue;
				}
				if (inverted != cmp < 0) {
					comparable = false;
				}
			}
			if (comparable && inverted) {
				swapLLLBasis(k, coefficients, basis, orthogonal, values, kmax, field, module);
				k = Math.max(1, k - 1);
			} else if (comparable) {
				for (int j = k - 2; j >= 0; j--) {
					reduceLLLBasis(k, j, coefficients, basis, field, module);
				}
				k++;
			} else {
				// We want to have kValue >= kValue or kValue < adjustedValue
				// by adjusting k with a unit.
				// unitlattice * x >= adjustedValue - kValue with x in Z^n
//				Vector<IntE> unitCoefficients = field.logRepresentationSpace().latticeLinearProgram(polytope, optimize,
//						integerLattice);
//				NFE unit = unitLattice.fromVector(unitCoefficients);
//				basis.set(k, module.scalarMultiply(unit, basis.get(k)));
			}
		}
		return basis;
	}

	private void reduceLLLBasis(int k, int j, List<List<NFE>> coefficients, List<Vector<NFE>> basis, NumberField field,
			FreeModule<NFE> module) {
		NFE m = field.minkowskiEmbeddingSpace()
				.closestLatticePoint(field.minkowskiEmbedding(coefficients.get(k).get(j)), field.maximalOrder());
		basis.set(k, module.subtract(basis.get(k), module.scalarMultiply(m, basis.get(j))));
		coefficients.get(k).set(j, field.subtract(coefficients.get(k).get(j), m));
		for (int i = 0; i < j; i++) {
			coefficients.get(k).set(i,
					field.subtract(coefficients.get(k).get(i), field.multiply(m, coefficients.get(j).get(i))));
		}
	}

	private void swapLLLBasis(int k, List<List<NFE>> coefficients, List<Vector<NFE>> basis,
			List<Vector<NFE>> orthogonal, List<NFE> values, int kmax, NumberField field, FreeModule<NFE> module) {
		Vector<NFE> tmp = basis.get(k);
		basis.set(k, basis.get(k - 1));
		basis.set(k - 1, tmp);
		for (int j = 0; j < k - 1; j++) {
			NFE c = coefficients.get(k).get(j);
			coefficients.get(k).set(j, coefficients.get(k - 1).get(j));
			coefficients.get(k - 1).set(j, c);
		}
		NFE c = coefficients.get(k).get(k - 1);
		NFE value = field.add(values.get(k), field.multiply(c, c, values.get(k - 1)));
		coefficients.get(k).set(k - 1, field.divide(field.multiply(c, values.get(k - 1)), value));
		values.set(k, field.divide(field.multiply(values.get(k - 1), values.get(k)), value));
		values.set(k - 1, value);
		Vector<NFE> tmpOrthogonal = orthogonal.get(k - 1);
		orthogonal.set(k - 1, module.add(orthogonal.get(k), module.scalarMultiply(c, orthogonal.get(k - 1))));
		orthogonal.set(k, module.subtract(tmpOrthogonal,
				module.scalarMultiply(coefficients.get(k).get(k - 1), orthogonal.get(k - 1))));
		for (int i = k + 1; i <= kmax; i++) {
			NFE m = coefficients.get(i).get(k);
			coefficients.get(i).set(k, field.subtract(coefficients.get(i).get(k - 1), field.multiply(m, c)));
			coefficients.get(i).set(k - 1,
					field.add(m, field.multiply(coefficients.get(k).get(k - 1), coefficients.get(i).get(k))));
		}
	}

}
