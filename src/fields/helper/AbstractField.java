package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.polynomials.Monomial;
import fields.polynomials.UnivariatePolynomialRing;

public abstract class AbstractField<T extends Element<T>> implements Field<T>, Ring<T> {
	private UnivariatePolynomialRing<T> ring = new UnivariatePolynomialRing<T>(this);

	@Override
	public Group<T> getAdditiveGroup() {
		return new Group<T>() {

			@Override
			public T getRandomElement() {
				return this.getRandomElement();
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return this.getNumberOfElements();
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractField.this.iterator();
			}

			@Override
			public T neutral() {
				return AbstractField.this.zero();
			}

			@Override
			public T inverse(T t) {
				return AbstractField.this.negative(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractField.this.add(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return this.isFinite();
			}

		};
	}

	@Override
	public Group<T> getMultiplicativeGroup() {
		return new Group<T>() {

			@Override
			public T getRandomElement() {
				T t;
				do {
					t = this.getRandomElement();
				} while (t.equals(AbstractField.this.zero()));
				return t;
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return this.getNumberOfElements().subtract(BigInteger.ONE);
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractField.this.getNonZeroElements().iterator();
			}

			@Override
			public T neutral() {
				return AbstractField.this.one();
			}

			@Override
			public T inverse(T t) {
				return AbstractField.this.inverse(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractField.this.multiply(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return AbstractField.this.isFinite();
			}
		};
	}

	@Override
	public T add(T t1, T t2, T t3) {
		return this.add(this.add(t1, t2), t3);
	}

	@Override
	public T subtract(T minuend, T subtrahend) {
		return this.add(minuend, this.negative(subtrahend));
	}

	@Override
	public T getInteger(int n) {
		return multiply(n, this.one());
	}

	@Override
	public T getInteger(BigInteger n) {
		return multiply(n, this.one());
	}

	@Override
	public T multiply(int n, T t) {
		return this.multiply(BigInteger.valueOf(n), t);
	}

	@Override
	public T multiply(BigInteger n, T t) {
		T result = this.zero();
		if (n.signum() < 0) {
			n = n.negate();
			t = this.negative(t);
		}
		for (int i = 0; i <= n.bitLength(); i++) {
			if (n.testBit(i))
				result = this.add(result, t);
			t = this.add(t, t);
		}
		return result;
	}

	@Override
	public T multiply(T t1, T t2, T t3) {
		return this.multiply(t1, this.multiply(t2, t3));
	}

	@Override
	public T multiply(int n, T t1, T t2) {
		return this.multiply(n, this.multiply(t1, t2));
	}

	@Override
	public T multiply(BigInteger n, T t1, T t2) {
		return this.multiply(n, this.multiply(t1, t2));
	}

	@Override
	public T multiply(int n, T t1, T t2, T t3) {
		return this.multiply(n, this.multiply(t1, t2, t3));
	}

	@Override
	public T multiply(BigInteger n, T t1, T t2, T t3) {
		return this.multiply(n, this.multiply(t1, t2, t3));
	}

	@Override
	public T divide(T dividend, T divisor) {
		return this.multiply(dividend, this.inverse(divisor));
	}

	@Override
	public T power(T t, int n) {
		return this.power(t, BigInteger.valueOf(n));
	}

	@Override
	public T power(T t, BigInteger n) {
		T result = this.one();
		if (n.signum() < 0) {
			n = n.negate();
			t = this.inverse(t);
		}
		for (int i = 0; i <= n.bitLength(); i++) {
			if (n.testBit(i))
				result = this.multiply(result, t);
			t = this.multiply(t, t);
		}
		return result;
	}

	public BigInteger discreteLogarithm(T base, T t) {
		if (!this.isFinite()) {
			throw new UnsupportedOperationException("not implemented");
		}
		throw new UnsupportedOperationException("not implemented");
	}

	public Iterable<T> getNonZeroElements() throws InfinityException {
		return new Iterable<T>() {
			private Iterable<T> elements = AbstractField.this;

			@Override
			public Iterator<T> iterator() {
				return new Iterator<T>() {
					private Iterator<T> it = elements.iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public T next() {
						T result = null;
						do {
							if (it.hasNext())
								result = it.next();
						} while (result.equals(AbstractField.this.zero()));
						return result;
					}

					@Override
					public void remove() {
						it.remove();
					}
				};
			}

		};
	}

	@Override
	public boolean isUnit(T t) {
		return !t.equals(this.zero());
	}
	
	@Override
	public T projectToUnit(T t) {
		if (t.equals(zero())) {
			return one();
		}
		return t;
	}
	
	@Override
	public T upToUnit(T t) {
		if (t.equals(zero())) {
			return zero();
		}
		return one();
	}

	@Override
	public List<Polynomial<T>> factorization(Polynomial<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean hasRoots(Polynomial<T> t) {
		List<Polynomial<T>> factors = factorization(t);
		for (Polynomial<T> factor : factors) {
			if (factor.degree() == 1) {
				return true;
			}
		}
		return false;
	}

	@Override
	public List<T> roots(Polynomial<T> t) {
		List<T> result = new ArrayList<>();
		Monomial constant = ring.getMonomial(new int[] { 0 });
		Monomial linear = ring.getMonomial(new int[] { 1 });
		List<Polynomial<T>> factors = factorization(t);
		for (Polynomial<T> factor : factors) {
			if (factor.degree() == 1) {
				result.add(divide(negative(factor.coefficient(constant)), factor.coefficient(linear)));
			}
		}
		return result;
	}

	@Override
	public boolean hasRoot(T t, int n) {
		if (t.equals(this.zero())) {
			return true;
		}
		Monomial m = ring.getMonomial(new int[] { n });
		Monomial c = ring.getMonomial(new int[] { 0 });
		SortedMap<Monomial, T> coeff = new TreeMap<>();
		coeff.put(m, this.one());
		coeff.put(c, this.negative(t));
		Polynomial<T> p = ring.getPolynomial(coeff);
		return this.hasRoots(p);
	}

	@Override
	public List<T> roots(T t, int n) {
		if (t.equals(this.zero())) {
			return Collections.singletonList(this.zero());
		}
		Monomial m = ring.getMonomial(new int[] { n });
		Monomial c = ring.getMonomial(new int[] { 0 });
		SortedMap<Monomial, T> coeff = new TreeMap<>();
		coeff.put(m, this.one());
		coeff.put(c, this.negative(t));
		Polynomial<T> p = ring.getPolynomial(coeff);
		List<T> roots = new ArrayList<>();
		roots.addAll(this.roots(p));
		return roots;
	}

	@Override
	public boolean hasSqrt(T t) {
		return hasRoot(t, 2);
	}

	@Override
	public List<T> sqrt(T t) {
		return roots(t, 2);
	}

	@Override
	public T primitiveRoot() {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(T t) {
		return t.equals(zero());
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public boolean isDivisible(T dividend, T divisor) {
		if (!divisor.equals(zero())) {
			return true;
		}
		return dividend.equals(zero());
	}

	@Override
	public List<T> quotientAndRemainder(T dividend, T divisor) {
		List<T> result = new ArrayList<T>();
		if (dividend.equals(zero()) && divisor.equals(zero())) {
			result.add(this.zero());
		} else {
			result.add(this.divide(dividend, divisor));
		}
		result.add(this.zero());
		return result;
	}
	
	@Override
	public T remainder(T dividend, T divisor) {
		return divisor.equals(zero()) ? dividend : zero();
	}

	@Override
	public BigInteger euclidMeasure(T t) {
		return BigInteger.ZERO;
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return this.getNumberOfElements().subtract(BigInteger.ONE);
	}

	@Override
	public Iterable<T> getUnits() {
		return this.getNonZeroElements();
	}

	public T gcd(T t1, T t2) {
		return this.extendedEuclidean(t1, t2).get(0);
	}

	public T lcm(T t1, T t2) {
		return divide(multiply(t1, t2), gcd(t1, t2));
	}

	public List<T> extendedEuclidean(T t1, T t2) {
		List<T> result = new ArrayList<T>();
		if (t1.equals(this.zero()) && t2.equals(this.zero())) {
			result.add(this.zero());
			result.add(this.zero());
			result.add(this.zero());
			return result;
		}
		result.add(this.one());
		if (t1.equals(this.zero())) {
			result.add(this.zero());
			result.add(this.inverse(t2));
			return result;
		}
		if (t2.equals(this.zero()) || this.characteristic().equals(BigInteger.TWO)) {
			result.add(this.inverse(t1));
			result.add(this.zero());
			return result;
		}
		result.add(this.inverse(this.multiply(2, t1)));
		result.add(this.inverse(this.multiply(2, t2)));
		return result;
	}

	@Override
	public UnivariatePolynomialRing<T> getUnivariatePolynomialRing() {
		return ring;
	}
	
	@Override
	public int krullDimension() {
		return 0;
	}

	public Ideal<T> getIdeal(List<T> generators) {
		return new FieldIdeal(generators);
	}

	public Ideal<T> add(Ideal<T> t1, Ideal<T> t2) {
		return getIdeal(Collections.singletonList(
				t1.generators().get(0).equals(zero()) && t2.generators().get(0).equals(zero()) ? zero() : one()));
	}

	public Ideal<T> multiply(Ideal<T> t1, Ideal<T> t2) {
		return getIdeal(Collections.singletonList(
				t1.generators().get(0).equals(zero()) || t2.generators().get(0).equals(zero()) ? zero() : one()));
	}
	
	public Ideal<T> intersect(Ideal<T> t1, Ideal<T> t2) {
		return getIdeal(Collections.singletonList(
				t1.generators().get(0).equals(zero()) || t2.generators().get(0).equals(zero()) ? zero() : one()));
	}

	public class FieldIdeal extends AbstractIdeal<T> {
		private boolean isZero;

		public FieldIdeal(List<T> generators) {
			super(AbstractField.this);
			isZero = true;
			for (T g : generators) {
				if (!g.equals(AbstractField.this.zero())) {
					isZero = false;
					break;
				}
			}
		}

		@Override
		public boolean isPrime() {
			return true;
		}

		@Override
		public boolean isMaximal() {
			return isZero;
		}

		@Override
		public List<T> generators() {
			return Collections.singletonList(isZero ? AbstractField.this.zero() : AbstractField.this.one());
		}

		@Override
		public List<T> generate(T t) {
			return Collections.singletonList(isZero ? zero() : t);
		}

		@Override
		public T residue(T t) {
			return isZero ? t : zero();
		}

		@Override
		public boolean isFinite() {
			return isZero || AbstractField.this.isFinite();
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			return isZero ? BigInteger.ONE : AbstractField.this.getNumberOfElements();
		}

		@Override
		public boolean contains(T t) {
			return !isZero || t.equals(zero());
		}
	}
}
