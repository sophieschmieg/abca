package fields.floatingpoint;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.ValueField;
import fields.polynomials.Monomial;

public class Complex extends AbstractField<ComplexNumber> implements ValueField<ComplexNumber> {
	private static Reals r = Reals.r();
	private static Complex c = new Complex();

	public static class ComplexNumber extends AbstractElement<ComplexNumber> {
		private Real realPart;
		private Real complexPart;

		public ComplexNumber(double realPart, double complexPart) {
			this.realPart = new Real(realPart);
			this.complexPart = new Real(complexPart);
		}

		public ComplexNumber(Real realPart, Real complexPart) {
			this.realPart = realPart;
			this.complexPart = complexPart;
		}

		@Override
		public int compareTo(ComplexNumber o) {
			int realComp = realPart.compareTo(o.realPart);
			if (realComp != 0) {
				return realComp;
			}
			return complexPart.compareTo(o.complexPart);
		}

		public Real realPart() {
			return realPart;
		}

		public Real complexPart() {
			return complexPart;
		}

		public String toString() {
			return realPart.toString() + " + " + complexPart.toString() + "i";
		}

	}

	private Complex() {
	}

	public static Complex c() {
		return c;
	}

	@Override
	public ComplexNumber zero() {
		return new ComplexNumber(r.zero(), r.zero());
	}

	@Override
	public ComplexNumber one() {
		return new ComplexNumber(r.one(), r.zero());
	}

	public ComplexNumber i() {
		return new ComplexNumber(r.zero(), r.one());
	}

	public double arg(ComplexNumber t) {
		return Math.atan2(t.complexPart.doubleValue(), t.realPart.doubleValue());
	}

	@Override
	public List<ComplexNumber> roots(ComplexNumber t, int n) {
		if (t.equals(zero())) {
			return Collections.nCopies(n, zero());
		}
		double arg = arg(t) / n;
		double value = Math.pow(value(t), 1.0 / n);
		List<ComplexNumber> roots = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			roots.add(new ComplexNumber(value * Math.cos(arg + 2 * i * Math.PI / n),
					value * Math.sin(arg + 2 * i * Math.PI / n)));
		}
		return roots;
	}

	@Override
	public ComplexNumber getInteger(BigInteger t) {
		return getEmbedding(r.getInteger(t));
	}

	public ComplexNumber getEmbedding(IntE t) {
		return getEmbedding(r.getEmbedding(t));
	}

	public ComplexNumber getEmbedding(Fraction t) {
		return getEmbedding(r.getEmbedding(t));
	}

	public ComplexNumber getEmbedding(Real t) {
		return new ComplexNumber(t, r.zero());
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public ComplexNumber add(ComplexNumber t1, ComplexNumber t2) {
		return new ComplexNumber(r.add(t1.realPart, t2.realPart), r.add(t1.complexPart, t2.complexPart));
	}

	@Override
	public ComplexNumber negative(ComplexNumber t) {
		return new ComplexNumber(r.negative(t.realPart), r.negative(t.complexPart));
	}

	@Override
	public ComplexNumber multiply(ComplexNumber t1, ComplexNumber t2) {
		return new ComplexNumber(
				t1.realPart.doubleValue() * t2.realPart.doubleValue()
						- t1.complexPart.doubleValue() * t2.complexPart.doubleValue(),
				t1.realPart.doubleValue() * t2.complexPart.doubleValue()
						+ t1.complexPart.doubleValue() * t2.realPart.doubleValue());
	}

	@Override
	public ComplexNumber inverse(ComplexNumber t) {
		ComplexNumber conjugate = conjugate(t);
		double denom = norm(t).doubleValue();
		return new ComplexNumber(conjugate.realPart.doubleValue() / denom, conjugate.complexPart.doubleValue() / denom);
	}

	@Override
	public ComplexNumber divide(ComplexNumber dividend, ComplexNumber divisor) {
		double denom = norm(divisor).doubleValue();
		ComplexNumber unnormed = multiply(dividend, conjugate(divisor));
		return new ComplexNumber(unnormed.realPart.doubleValue() / denom, unnormed.complexPart.doubleValue() / denom);
	}

	public ComplexNumber conjugate(ComplexNumber t) {
		return new ComplexNumber(t.realPart, r.negative(t.complexPart));
	}

	@Override
	public ComplexNumber getRandomElement() {
		return new ComplexNumber(r.getRandomElement(), r.getRandomElement());
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
	public Iterator<ComplexNumber> iterator() {
		throw new InfinityException();
	}

	public Real norm(ComplexNumber t) {
		return multiply(t, conjugate(t)).realPart;
	}

	@Override
	public double value(ComplexNumber t) {
		return Math.sqrt(norm(t).doubleValue());
	}

	public ComplexNumber exp(ComplexNumber t) {
		return multiply(getEmbedding(new Real(Math.exp(t.realPart.doubleValue()))),
				new ComplexNumber(Math.cos(t.complexPart.doubleValue()), Math.sin(t.complexPart.doubleValue())));
	}

	public ComplexNumber log(ComplexNumber t) {
		return new ComplexNumber(Math.log(value(t)), arg(t));
	}

	public ComplexNumber findRoot(Polynomial<ComplexNumber> t) {
		PolynomialRing<ComplexNumber> ring = getUnivariatePolynomialRing();
		Polynomial<ComplexNumber> derivative = ring.derivative(t, 1);
		double scale = 0.0;
		double lcScale = value(t.leadingCoefficient());
		for (Monomial m : t.monomials()) {
			if (m.degree() == t.degree()) {
				continue;
			}
			scale += Math.pow(value(t.coefficient(m)) / lcScale, 1.0 / (t.degree() - m.degree()));
		}
		int maxNumIterations = 1000;
		while (true) {
			ComplexNumber x = multiply(getEmbedding(new Real(scale)), getRandomElement());
			for (int i = 0; i < maxNumIterations; i++) {
				x = subtract(x, divide(ring.evaluate(t, x), ring.evaluate(derivative, x)));
				if (norm(ring.evaluate(t, x)).doubleValue() < Reals.EPSILON * Reals.EPSILON) {
					return x;
				}
			}
		}
	}

	public List<Polynomial<ComplexNumber>> factorization(Polynomial<ComplexNumber> t) {
		List<Polynomial<ComplexNumber>> squareFree = getUnivariatePolynomialRing().squareFreeFactorization(t);
		List<Polynomial<ComplexNumber>> result = new ArrayList<>();
		for (Polynomial<ComplexNumber> squareFreeFactor : squareFree) {
			result.addAll(squareFreeFactorization(squareFreeFactor));
		}
		return result;
	}

	public List<Polynomial<ComplexNumber>> squareFreeFactorization(Polynomial<ComplexNumber> t) {
		Polynomial<ComplexNumber> f = t;
		List<Polynomial<ComplexNumber>> factors = new ArrayList<>();
		if (t.degree() <= 0) {
			return Collections.singletonList(t);
		}
		PolynomialRing<ComplexNumber> ring = getUnivariatePolynomialRing();
		while (f.degree() > 0) {
			ComplexNumber root = findRoot(f);
			Polynomial<ComplexNumber> factor = ring.subtract(ring.getVar(1), ring.getEmbedding(root));
			factors.add(factor);
			f = ring.quotientAndRemainder(f, factor).get(0);
		}
		return factors;
	}
}
