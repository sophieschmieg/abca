package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.NTTRing.NTT;
import fields.helper.ProductRing.ProductElement;
import fields.helper.ProductRing.ProductIdeal;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.Vector;

public class NTTRing<T extends Element<T>, R extends Ring<T>> extends AbstractRing<NTT<T>> {
	public static class NTT<T extends Element<T>> extends AbstractElement<NTT<T>> {
		private ProductElement<T> asProduct;

		private NTT(ProductElement<T> asProduct) {
			this.asProduct = asProduct;
		}

		private ProductElement<T> asProduct() {
			return asProduct;
		}

		@Override
		public String toString() {
			return asProduct.toString();
		}

		@Override
		public int compareTo(NTT<T> o) {
			return asProduct.compareTo(o.asProduct);
		}

	}

	private R base;
	private int degree;
	private ProductRing<T, R> product;
	private T inverseDegree;
	private List<T> rootPowers;
	private List<T> inverseRootPowers;

	public NTTRing(R base, int degree, T rootOfUnity) {
		this.degree = degree;
		this.base = base;
		if (!base.power(rootOfUnity, degree).equals(base.getInteger(-1))) {
			throw new ArithmeticException(rootOfUnity + " is not a " + (2 * degree) + " root of unity!");
		}
		T inverseRootOfUnity = base.inverse(rootOfUnity);
		this.inverseDegree = base.inverse(base.getInteger(degree));
		this.product = ProductRing.power(base, degree);
		this.rootPowers = new ArrayList<>();
		this.inverseRootPowers = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			int exponent = bitReverse(i, degree);
			rootPowers.add(base.power(rootOfUnity, exponent));
			inverseRootPowers.add(base.power(inverseRootOfUnity, exponent));
		}
	}

	@Override
	public String toString() {
		return base + "[X]/(X^" + degree + " + 1)";
	}

	private NTT<T> getElement(ProductElement<T> t) {
		return new NTT<>(t);
	}

	@Override
	public NTT<T> zero() {
		return getElement(product.zero());
	}

	@Override
	public NTT<T> one() {
		return getElement(product.one());
	}

	private int bitReverse(int index, int degree) {
		int inputPower = degree >> 1;
		int outputPower = 1;
		int result = 0;
		while (inputPower > 0) {
			result |= ((index & inputPower) != 0) ? outputPower : 0;
			inputPower >>= 1;
			outputPower <<= 1;
		}
		return result;
	}

	public NTT<T> fromPolynomial(UnivariatePolynomial<T> t) {
		List<T> coeffs = new ArrayList<>();
		coeffs.addAll(t.coefficients());
		while (coeffs.size() < degree) {
			coeffs.add(base.zero());
		}
		while (coeffs.size() > degree) {
			int index = coeffs.size() - 1;
			int newIndex = index % degree;
			int power = (index - newIndex) / degree;
			T last = base.multiply(power % 2 == 0 ? base.one() : base.getInteger(-1), coeffs.remove(index));
			coeffs.set(newIndex, base.add(coeffs.get(newIndex), last));
		}
		int offset = degree;
		for (int step = 1; step < degree; step *= 2) {
			offset /= 2;
			int k = 0;
			for (int i = 0; i < step; i++) {
				T stepRoot = rootPowers.get(i + step);
				for (int j = k; j < k + offset; j++) {
					T odd = base.multiply(stepRoot, coeffs.get(j + offset));
					T even = coeffs.get(j);
					coeffs.set(j, base.add(odd, even));
					coeffs.set(j + offset, base.subtract(even, odd));
				}
				k += 2 * offset;
			}
		}
		return getElement(product.getElement(coeffs));
	}

	public UnivariatePolynomial<T> asPolynomial(NTT<T> t) {
		List<T> coeffs = new ArrayList<>();
		coeffs.addAll(t.asProduct().values());
		int offset = degree;
		for (int step = 1; step < degree; step *= 2) {
			offset /= 2;
			int k = 0;
			for (int i = 0; i < offset; i++) {
				T stepRoot = inverseRootPowers.get(i + offset);
				for (int j = k; j < k + step; j++) {
					T odd = coeffs.get(j + step);
					T even = coeffs.get(j);
					coeffs.set(j, base.add(odd, even));
					coeffs.set(j + step, base.multiply(stepRoot, base.subtract(even, odd)));
				}
				k += 2 * step;
			}
		}
		UnivariatePolynomialRing<T> polynomialRing = base.getUnivariatePolynomialRing();
		return polynomialRing.multiply(inverseDegree, polynomialRing.getPolynomial(coeffs));
	}

	public NTT<T> fromList(List<T> t) {
		return getElement(product.getElement(t));
	}

	public List<T> asList(NTT<T> t) {
		return t.asProduct().values();
	}

	@Override
	public BigInteger characteristic() {
		return base.characteristic();
	}

	@Override
	public NTT<T> add(NTT<T> t1, NTT<T> t2) {
		return getElement(product.add(t1.asProduct(), t2.asProduct()));
	}

	@Override
	public NTT<T> negative(NTT<T> t) {
		return getElement(product.negative(t.asProduct()));
	}

	@Override
	public NTT<T> multiply(NTT<T> t1, NTT<T> t2) {
		return getElement(product.multiply(t1.asProduct(), t2.asProduct()));
	}

	@Override
	public boolean isUnit(NTT<T> t) {
		return product.isUnit(t.asProduct());
	}

	@Override
	public NTT<T> inverse(NTT<T> t) {
		return getElement(product.inverse(t.asProduct()));
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isIntegral() {
		return product.isIntegral();
	}

	@Override
	public boolean isReduced() {
		return product.isReduced();
	}

	@Override
	public boolean isIrreducible() {
		return product.isIrreducible();
	}

	@Override
	public boolean isZeroDivisor(NTT<T> t) {
		return product.isZeroDivisor(t.asProduct());
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<NTT<T>, NTT<T>> uniqueFactorization(NTT<T> t) {
		throw new ArithmeticException("Not a UFD!");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return true;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(NTT<T> dividend, NTT<T> divisor) {
		return product.isDivisible(dividend.asProduct(), divisor.asProduct());
	}

	@Override
	public QuotientAndRemainderResult<NTT<T>> quotientAndRemainder(NTT<T> dividend, NTT<T> divisor) {
		QuotientAndRemainderResult<ProductElement<T>> qr = product.quotientAndRemainder(dividend.asProduct(),
				divisor.asProduct());
		return new QuotientAndRemainderResult<>(getElement(qr.getQuotient()), getElement(qr.getRemainder()));
	}

	@Override
	public BigInteger euclidMeasure(NTT<T> t) {
		throw new ArithmeticException("not a euclidean ring!");
	}

	@Override
	public NTT<T> projectToUnit(NTT<T> t) {
		return getElement(product.projectToUnit(t.asProduct()));
	}

	@Override
	public Iterable<NTT<T>> getUnits() {
		return new Iterable<>() {

			@Override
			public Iterator<NTT<T>> iterator() {
				return new Iterator<>() {
					private Iterator<ProductElement<T>> it = product.getUnits().iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public NTT<T> next() {
						return getElement(it.next());
					}

				};
			}
		};
	}

	@Override
	public int krullDimension() {
		return 0;
	}

	@Override
	public IdealResult<NTT<T>, NTTIdeal> getIdealWithTransforms(List<NTT<T>> generators) {
		List<ProductElement<T>> productGenerators = new ArrayList<>();
		for (NTT<T> generator : generators) {
			productGenerators.add(generator.asProduct());
		}
		IdealResult<ProductElement<T>, ProductIdeal<T, R>> ideal = product.getIdealWithTransforms(productGenerators);
		List<List<NTT<T>>> expressions = new ArrayList<>();
		List<Vector<NTT<T>>> syzygies = new ArrayList<>();
		for (List<ProductElement<T>> expression : ideal.getGeneratorExpressions()) {
			List<NTT<T>> nttExpression = new ArrayList<>();
			for (ProductElement<T> t : expression) {
				nttExpression.add(getElement(t));
			}
			expressions.add(nttExpression);
		}
		for (Vector<ProductElement<T>> syzygy : ideal.getSyzygies()) {
			List<NTT<T>> nttSyzygy = new ArrayList<>();
			for (ProductElement<T> t : syzygy.asList()) {
				nttSyzygy.add(getElement(t));
			}
			syzygies.add(new Vector<>(nttSyzygy));
		}
		return new IdealResult<>(expressions, generators, new NTTIdeal(ideal.getIdeal()), syzygies);
	}

	@Override
	public NTTIdeal intersect(Ideal<NTT<T>> t1, Ideal<NTT<T>> t2) {
		return new NTTIdeal(product.intersect(((NTTIdeal) t1).asProductIdeal, ((NTTIdeal) t2).asProductIdeal));
	}

	@Override
	public NTTIdeal radical(Ideal<NTT<T>> t) {
		return new NTTIdeal(product.radical(((NTTIdeal) t).asProductIdeal));
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public NTT<T> getRandomElement() {
		return getElement(product.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return product.getNumberOfElements();
	}

	@Override
	public Iterator<NTT<T>> iterator() {
		return new Iterator<>() {
			private Iterator<ProductElement<T>> it = product.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public NTT<T> next() {
				return getElement(it.next());
			}
		};
	}

	public class NTTIdeal extends AbstractIdeal<NTT<T>> {
		private ProductIdeal<T, R> asProductIdeal;

		private NTTIdeal(ProductIdeal<T, R> asProductIdeal) {
			super(NTTRing.this);
			this.asProductIdeal = asProductIdeal;
		}

		@Override
		public boolean isPrimary() {
			return asProductIdeal.isPrimary();
		}

		@Override
		public boolean isPrime() {
			return asProductIdeal.isPrime();
		}

		@Override
		public boolean isMaximal() {
			return asProductIdeal.isMaximal();
		}

		@Override
		public List<NTT<T>> generators() {
			List<NTT<T>> generators = new ArrayList<>();
			for (ProductElement<T> generator : asProductIdeal.generators()) {
				generators.add(getElement(generator));
			}
			return generators;
		}

		@Override
		public List<NTT<T>> generate(NTT<T> t) {
			List<NTT<T>> result = new ArrayList<>();
			for (ProductElement<T> factor : asProductIdeal.generate(t.asProduct())) {
				result.add(getElement(factor));
			}
			return result;
		}

		@Override
		public NTT<T> residue(NTT<T> t) {
			return getElement(asProductIdeal.residue(t.asProduct()));
		}

		@Override
		public boolean contains(NTT<T> t) {
			return asProductIdeal.contains(t.asProduct());
		}

		@Override
		public boolean isFinite() {
			return true;
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			return asProductIdeal.getNumberOfElements();
		}
	}
}
