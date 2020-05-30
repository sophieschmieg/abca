package fields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.CoordinateRing.CoordinateRingElement;

public class CoordinateRing<T extends Element> extends AbstractAlgebra<Polynomial<T>, CoordinateRingElement<T>>
		implements Ring<CoordinateRingElement<T>> {
	private PolynomialRing<T> ring;
	private PolynomialRing<T>.Ideal ideal;

	public CoordinateRing(PolynomialRing<T> ring, PolynomialRing<T>.Ideal ideal) {
		this.ring = ring;
		this.ideal = ideal;
	}

	@Override
	public String toString() {
		return this.ring.toString() + "/" + this.ideal.toString();
	}

	public PolynomialRing<T> getPolynomialRing() {
		return ring;
	}
	
	public Ring<Polynomial<T>> getRing() {
		return ring;
	}
	
	public PolynomialRing<T>.Ideal getIdeal() {
		return ideal;
	}
	
	public boolean isFree() {
		throw new ArithmeticException();
	}
	
	public List<CoordinateRingElement<T>> getGenerators() {
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		for (int i = 1; i <= this.ring.getNumVars(); i++) {
			result.add(getVar(i));
		}
		return result;
	}
	
	public CoordinateRingElement<T> getVar(int var) {
		return this.getEmbedding(this.ring.getVar(var));
	}
	
	public CoordinateRingElement<T> getEmbedding(Polynomial<T> t) {
		return new CoordinateRingElement<T>(this.ideal, t);
	}

	private CoordinateRingElement<T> getEmbedding(Polynomial<T> t, boolean reduce) {
		return new CoordinateRingElement<T>(this.ideal, t, reduce);
	}

	public CoordinateRingElement<T> getEmbedding(T t) {
		return this.getEmbedding(this.ring.getEmbedding(t), false);
	}

	@Override
	public CoordinateRingElement<T> zero() {
		return this.getEmbedding(this.ring.zero(), false);
	}

	@Override
	public CoordinateRingElement<T> one() {
		return this.getEmbedding(this.ring.one(), false);
	}

	@Override
	public CoordinateRingElement<T> add(CoordinateRingElement<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.add(t1.polynomial, t2.polynomial), false);
	}

	@Override
	public CoordinateRingElement<T> negative(CoordinateRingElement<T> t) {
		return this.getEmbedding(this.ring.negative(t.polynomial), false);
	}

	@Override
	public CoordinateRingElement<T> multiply(CoordinateRingElement<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1.polynomial, t2.polynomial));
	}
	
	public CoordinateRingElement<T> multiply(T t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}
	
	public CoordinateRingElement<T> scalarMultiply(Polynomial<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	@Override
	public CoordinateRingElement<T> inverse(CoordinateRingElement<T> t) {
		PolynomialRing<T> r = ring.addVariableWithElimination(1);
		Polynomial<T> s = r.getVar(1);
		Polynomial<T> f = r.getEmbeddingShift(t.getElement(), 1);
		Polynomial<T> inv = r.subtract(r.one(), r.multiply(s, f));
		List<Polynomial<T>> basis = new ArrayList<Polynomial<T>>();
		basis.add(inv);
		for (Polynomial<T> b : ideal.getBasis()) {
			basis.add(r.getEmbeddingShift(b, 1));
		}
		PolynomialRing<T>.Ideal i = r.getIdeal(basis);
		Polynomial<T> inverted = i.getBasis().first();
		int[] expectedExponents = new int[r.getNumVars()];
		expectedExponents[0] = 1;
		if (!inverted.getLeadingTerm().equals(new Polynomial.Monomial(expectedExponents))) {
			throw new ArithmeticException("No inverse found!");
		}
		Polynomial<T> shifted = ring.getEmbeddingShift(r.add(s, r.negative(inverted)), -1);
		return this.getEmbedding(shifted);
	}

	@Override
	public boolean isUnit(CoordinateRingElement<T> t) {
		PolynomialRing<T>.Ideal tIdeal = this.ring.getIdeal(Collections.singletonList(t.polynomial));
		PolynomialRing<T>.Ideal degenerate = this.ring.getIdeal(Collections.singletonList(this.ring.one()));
		return this.ideal.add(tIdeal).equalsIdeal(degenerate);
	}

	@Override
	public CoordinateRingElement<T> getRandomElement() {
		throw new UnsupportedOperationException();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return BigInteger.valueOf(-1);
	}

	@Override
	public Iterator<CoordinateRingElement<T>> iterator() {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isFinite() {
		return this.ring.getField().isFinite() && this.ideal.dimension() == 0;
	}

	@Override
	public BigInteger characteristic() {
		return this.ring.characteristic();
	}

	@Override
	public boolean isIntegral() {
		return false; // TODO check whether integral
	}
	
	@Override
	public boolean isZeroDivisor(CoordinateRingElement<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}
	
	@Override
	public boolean isDivisible(CoordinateRingElement<T> dividend,
			CoordinateRingElement<T> divisor) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<CoordinateRingElement<T>> quotientAndRemainder(CoordinateRingElement<T> dividend,
			CoordinateRingElement<T> divisor) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public BigInteger euclidMeasure(CoordinateRingElement<T> t) {
		return null;
	}

	@Override
	public Iterable<CoordinateRingElement<T>> getUnits() {
		throw new UnsupportedOperationException();
	}

	public CoordinateRingElement<T> substitute(CoordinateRingElement<T> t, List<CoordinateRingElement<T>> values) {
		CoordinateRingElement<T> result = this.zero();
		for (Polynomial.Monomial m : t.polynomial.getCoefficients().keySet()) {
			CoordinateRingElement<T> value = this.getEmbedding(t.polynomial.getCoefficient(m));
			int exponents[] = m.getExponents();
			for (int i = 0; i < exponents.length; i++) {
				value = this.multiply(value, this.power(values.get(i), exponents[i]));
			}
			result = this.add(result, value);
		}
		return result;
	}

	public static class CoordinateRingElement<T extends Element> implements Element {
		private Polynomial<T> polynomial;

		public CoordinateRingElement(PolynomialRing<T>.Ideal ideal, Polynomial<T> polynomial) {
			this.polynomial = ideal.residue(polynomial);
		}

		private CoordinateRingElement(PolynomialRing<T>.Ideal ideal, Polynomial<T> polynomial, boolean reduce) {
			if (reduce) {
				this.polynomial = ideal.residue(polynomial);
			} else {
				this.polynomial = polynomial;
			}
		}

		public boolean equals(Object o) {
			if (!(o instanceof CoordinateRing.CoordinateRingElement))
				return false;
			@SuppressWarnings("unchecked")
			CoordinateRingElement<T> t = (CoordinateRingElement<T>) o;
			return this.polynomial.equals(t.polynomial);
		}

		public String toString() {
			return this.polynomial.toString();
		}

		public Polynomial<T> getElement() {
			return this.polynomial;
		}

		@Override
		public int compareTo(Element e) {
			@SuppressWarnings("unchecked")
			CoordinateRingElement<T> o = (CoordinateRingElement<T>) e;
			return this.polynomial.compareTo(o.polynomial);
		}
	}
}
