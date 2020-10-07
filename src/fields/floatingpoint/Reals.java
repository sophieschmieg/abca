package fields.floatingpoint;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.ValueField;
import util.MiscAlgorithms;

public class Reals extends AbstractField<Real> implements ValueField<Real> {
	private static Reals r = new Reals();
	public static double EPSILON = Math.pow(2.0, -50);

	public static class Real extends AbstractElement<Real> {
		private double value;

		public Real(double value) {
			this.value = value;
		}

		@Override
		public int compareTo(Real o) {
			if (Math.abs(value - o.value) < EPSILON) {
				return 0;
			}
			return value < o.value ? -1 : 1;
		}

		public double doubleValue() {
			return value;
		}

		public String toString() {
			return Double.toString(value);
		}
	}

	private Reals() {
	}

	public static Reals r() {
		return r;
	}

	@Override
	public Real zero() {
		return new Real(0.0);
	}

	@Override
	public Real one() {
		return new Real(1.0);
	}

	@Override
	public Real getInteger(BigInteger t) {
		return new Real(t.doubleValue());
	}

	public Real getEmbedding(IntE t) {
		return getInteger(t.getValue());
	}

	public Real getEmbedding(Fraction t) {
		return divide(getEmbedding(t.getNumerator()), getEmbedding(t.getDenominator()));
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}
	
	@Override
	public List<Real> roots(Real t, int n) {
		if (t.value == 0.0) {
			return Collections.nCopies(n, zero());
		}
		if (n % 2 == 0 && t.value < 0.0) {
			return Collections.emptyList();
		}
		Real root = new Real(Math.signum(t.value) * Math.pow(Math.abs(t.value), 1.0/n));
		List<Real> roots = new ArrayList<>();
		roots.add(root);
		if (n % 2 == 0) {
			roots.add(negative(root));
		}
		return roots;
	}

	@Override
	public Real add(Real t1, Real t2) {
		return new Real(t1.value + t2.value);
	}

	@Override
	public Real negative(Real t) {
		return new Real(-t.value);
	}

	@Override
	public Real multiply(Real t1, Real t2) {
		return new Real(t1.value * t2.value);
	}

	@Override
	public Real inverse(Real t) {
		return new Real(1.0 / t.value);
	}
	
	@Override
	public Real divide(Real dividend, Real divisor) {
		return new Real(dividend.value / divisor.value);
	}

	@Override
	public Real getRandomElement() {
		return new Real(new Random().nextGaussian());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<Real> iterator() {
		throw new InfinityException();
	}

	@Override
	public double value(Real t) {
		return Math.abs(t.value);
	}

	@Override
	public List<Polynomial<Real>> factorization(Polynomial<Real> t) {
		List<ComplexNumber> complexRoots = Complex.c().roots(MiscAlgorithms.mapPolynomial(t, new MathMap<>() {

			@Override
			public ComplexNumber evaluate(Real n) {
				return Complex.c().getEmbedding(n);
			}
		}, Complex.c().getUnivariatePolynomialRing()));
		List<Polynomial<Real>> factors = new ArrayList<>();
		PolynomialRing<Real> ring = getUnivariatePolynomialRing();
		for (ComplexNumber root : complexRoots) {
			if (root.complexPart().equals(zero())) {
				factors.add(ring.subtract(ring.getVar(1), ring.getEmbedding(root.realPart())));
			} else if (root.complexPart().doubleValue() > 0.0) {
				double value = Complex.c().value(root);
				double norm = value * value;
				factors.add(
						ring.add(ring.getVarPower(1, 2), ring.getEmbedding(r.multiply(-2, root.realPart()), new int[] { 1 }),
								ring.getEmbedding(new Real(norm))));
			}
		}
		return factors;
	}
}
