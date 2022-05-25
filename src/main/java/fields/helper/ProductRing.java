package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.ProductRing.ProductElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.Vector;

public class ProductRing<T extends Element<T>, R extends Ring<T>> extends AbstractRing<ProductElement<T>>
		implements Ring<ProductElement<T>> {
	public static class ProductElement<T extends Element<T>> extends AbstractElement<ProductElement<T>> {
		private List<T> values;

		private ProductElement(List<T> values) {
			this.values = values;
		}

		public List<T> values() {
			return values;
		}

		public T get(int index) {
			return values.get(index);
		}

		@Override
		public int compareTo(ProductElement<T> o) {
			if (values.size() != o.values.size()) {
				return values.size() - o.values.size();
			}
			for (int i = 0; i < values.size(); i++) {
				int cmp = get(i).compareTo(o.get(i));
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}
	}

	public static <T extends Element<T>, R extends Ring<T>> ProductRing<T, R> power(R ring, int power) {
		List<R> rings = new ArrayList<>();
		for (int i = 0; i < power; i++) {
			rings.add(ring);
		}
		return new ProductRing<>(rings);
	}

	private List<R> rings;

	public ProductRing(List<R> rings) {
		this.rings = rings;
	}

	public List<R> getRings() {
		return rings;
	}

	public R getRing(int index) {
		return rings.get(index);
	}

	public int numberOfFactors() {
		return rings.size();
	}

	public ProductElement<T> getElement(List<T> values) {
		return new ProductElement<>(values);
	}

	public ProductElement<T> inject(T element, int index) {
		List<T> values = new ArrayList<>();
		for (int i = 0; i < index; i++) {
			values.add(rings.get(i).zero());
		}
		values.add(element);
		for (int i = index + 1; i < numberOfFactors(); i++) {
			values.add(rings.get(i).zero());
		}
		return getElement(values);
	}

	public MathMap<T, ProductElement<T>> injection(int index) {
		return new MathMap<>() {
			@Override
			public ProductElement<T> evaluate(T t) {
				return inject(t, index);
			}
		};
	}

	public T project(ProductElement<T> t, int index) {
		return t.get(index);
	}

	public MathMap<ProductElement<T>, T> projection(int index) {
		return new MathMap<>() {
			@Override
			public T evaluate(ProductElement<T> t) {
				return project(t, index);
			}
		};
	}

	public ProductElement<T> constantInteger(IntE t) {
		List<T> values = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			values.add(rings.get(i).getInteger(t));
		}
		return getElement(values);
	}

	@Override
	public Exactness exactness() {
		if (numberOfFactors() == 0) {
			return Exactness.EXACT;
		}
		return rings.get(0).exactness();
	}

	@Override
	public ProductElement<T> getRandomElement() {
		List<T> values = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			values.add(rings.get(i).getRandomElement());
		}
		return getElement(values);
	}

	@Override
	public boolean isFinite() {
		for (R ring : rings) {
			if (!ring.isFinite()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		BigInteger number = BigInteger.ONE;
		for (R ring : rings) {
			number = number.multiply(ring.getNumberOfElements());
		}
		return number;
	}

	private Iterator<ProductElement<T>> iterator(List<? extends Iterable<T>> iterables) {
		if (iterables.isEmpty()) {
			return new Iterator<ProductElement<T>>() {
				private boolean used = false;

				@Override
				public boolean hasNext() {
					return !used;
				}

				@Override
				public ProductElement<T> next() {
					if (used) {
						throw new RuntimeException("Iterator at end!");
					}
					used = true;
					return zero();
				}
			};
		}
		return new Iterator<ProductElement<T>>() {
			private int limitIndex;
			private long limit;
			private long[] index;
			private List<T> current = new ArrayList<>();
			private List<Iterator<T>> its = null;

			@Override
			public boolean hasNext() {
				return current != null;
			}

			@Override
			public ProductElement<T> next() {
				if (its == null) {
					its = new ArrayList<>();
					limit = 0;
					index = new long[iterables.size()];
					limitIndex = iterables.size() - 1;
					for (int i = 0; i < iterables.size(); i++) {
						its.add(iterables.get(i).iterator());
						current.add(its.get(i).next());
						index[i] = 0;
					}
				}
				if (current == null) {
					throw new RuntimeException("Iterator at end!");
				}
				List<T> copy = new ArrayList<>();
				copy.addAll(current);
				ProductElement<T> next = getElement(copy);
				boolean broken = false;
				boolean hasNext = false;
				for (int i = 0; i < iterables.size(); i++) {
					if (its.get(i).hasNext()) {
						hasNext = true;
						break;
					}
				}
				while (!broken && hasNext) {
					for (int i = 0; i < iterables.size(); i++) {
						if (!its.get(i).hasNext()) {
							continue;
						}
						if (i == limitIndex) {
							continue;
						}
						if (i < limitIndex && index[i] >= limit) {
							continue;
						}
						if (i > limitIndex && index[i] >= limit - 1) {
							continue;
						}
						current.set(i, its.get(i).next());
						index[i]++;
						for (int j = 0; j < i; j++) {
							if (j == limitIndex) {
								continue;
							}
							its.set(j, iterables.get(j).iterator());
							current.set(j, its.get(j).next());
							index[j] = 0;
						}
						broken = true;
						break;
					}
					if (!broken) {
						if (limitIndex < iterables.size()) {
							limitIndex++;
							current.set(limitIndex, its.get(limitIndex).next());
							index[limitIndex]++;
							for (int i = 0; i < iterables.size(); i++) {
								if (i == limitIndex) {
									continue;
								}
								its.set(i, iterables.get(i).iterator());
								current.set(i, its.get(i).next());
								index[i] = 0;

							}
						} else {
							limit++;
							limitIndex = 0;
							current.set(0, its.get(0).next());
							index[0]++;
							for (int i = 1; i < iterables.size(); i++) {
								its.set(i, iterables.get(i).iterator());
								current.set(i, its.get(i).next());
								index[i] = 0;
							}
						}
						break;
					}
				}
				if (!hasNext) {
					current = null;
				}
				return next;
			}
		};
	}

	@Override
	public Iterator<ProductElement<T>> iterator() {
		return iterator(rings);
	}

	@Override
	public ProductElement<T> zero() {
		return constantInteger(Integers.z().zero());
	}

	@Override
	public ProductElement<T> one() {
		return constantInteger(Integers.z().one());
	}

	@Override
	public BigInteger characteristic() {
		Integers z = Integers.z();
		IntE characteristic = z.one();
		for (R ring : rings) {
			characteristic = z.lcm(characteristic, z.getInteger(ring.characteristic()));
		}
		return characteristic.getValue();
	}

	@Override
	public ProductElement<T> add(ProductElement<T> t1, ProductElement<T> t2) {
		List<T> values = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			values.add(rings.get(i).add(t1.get(i), t2.get(i)));
		}
		return getElement(values);
	}

	@Override
	public ProductElement<T> negative(ProductElement<T> t) {
		List<T> values = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			values.add(rings.get(i).negative(t.get(i)));
		}
		return getElement(values);
	}

	@Override
	public ProductElement<T> multiply(ProductElement<T> t1, ProductElement<T> t2) {
		List<T> values = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			values.add(rings.get(i).multiply(t1.get(i), t2.get(i)));
		}
		return getElement(values);
	}

	@Override
	public boolean isUnit(ProductElement<T> t) {
		for (int i = 0; i < numberOfFactors(); i++) {
			if (!rings.get(i).isUnit(t.get(i))) {
				return false;
			}
		}
		return true;
	}

	@Override
	public ProductElement<T> inverse(ProductElement<T> t) {
		List<T> values = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			values.add(rings.get(i).inverse(t.get(i)));
		}
		return getElement(values);
	}

	@Override
	public boolean isCommutative() {
		for (R ring : rings) {
			if (!ring.isCommutative()) {
				return false;
			}
		}
		return true;
	}

	private Optional<Integer> uniqueNonZeroRingIndex() {
		int result = -1;
		for (int i = 0; i < numberOfFactors(); i++) {
			R ring = rings.get(i);
			if (!ring.getNumberOfElements().equals(BigInteger.ONE)) {
				if (result != -1) {
					return Optional.empty();
				}
				result = i;
			}
		}
		return result == -1 ? Optional.empty() : Optional.of(result);

	}

	private Optional<R> uniqueNonZeroRing() {
		Optional<Integer> nonZero = uniqueNonZeroRingIndex();
		if (nonZero.isEmpty()) {
			return Optional.empty();
		}
		return Optional.of(rings.get(nonZero.get()));
	}

	@Override
	public boolean isIntegral() {
		Optional<R> nonZero = uniqueNonZeroRing();
		return nonZero.isPresent() && nonZero.get().isIntegral();
	}

	@Override
	public boolean isReduced() {
		for (R ring : rings) {
			if (!ring.isReduced()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isIrreducible() {
		Optional<R> nonZero = uniqueNonZeroRing();
		return nonZero.isPresent() && nonZero.get().isIrreducible();
	}

	@Override
	public boolean isZeroDivisor(ProductElement<T> t) {
		for (int i = 0; i < numberOfFactors(); i++) {
			if (rings.get(i).isZeroDivisor(t.get(i))) {
				return true;
			}
		}
		return false;
	}

	@Override
	public boolean isEuclidean() {
		Optional<R> nonZero = uniqueNonZeroRing();
		return nonZero.isPresent() && nonZero.get().isEuclidean();
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		Optional<R> nonZero = uniqueNonZeroRing();
		return nonZero.isPresent() && nonZero.get().isUniqueFactorizationDomain();
	}

	@Override
	public FactorizationResult<ProductElement<T>, ProductElement<T>> uniqueFactorization(ProductElement<T> t) {
		Optional<Integer> nonZero = uniqueNonZeroRingIndex();
		if (nonZero.isEmpty() || !rings.get(nonZero.get()).isUniqueFactorizationDomain()) {
			throw new ArithmeticException("Not a UFD!");
		}
		R ring = rings.get(nonZero.get());
		FactorizationResult<T, T> factorization = ring.uniqueFactorization(project(t, nonZero.get()));
		SortedMap<ProductElement<T>, Integer> result = new TreeMap<>();
		for (T prime : factorization.primeFactors()) {
			result.put(inject(prime, nonZero.get()), factorization.multiplicity(prime));
		}
		return new FactorizationResult<>(inject(factorization.getUnit(), nonZero.get()), result);
	}

	@Override
	public FactorizationResult<Polynomial<ProductElement<T>>, ProductElement<T>> factorization(
			UnivariatePolynomial<ProductElement<T>> t) {
		Optional<Integer> nonZero = uniqueNonZeroRingIndex();
		if (nonZero.isEmpty() || !rings.get(nonZero.get()).isUniqueFactorizationDomain()) {
			throw new ArithmeticException("Not a UFD!");
		}
		R ring = rings.get(nonZero.get());
		UnivariatePolynomialRing<T> polynomialRing = ring.getUnivariatePolynomialRing();
		FactorizationResult<Polynomial<T>, T> factorization = ring
				.factorization(polynomialRing.getEmbedding(t, projection(nonZero.get())));
		SortedMap<Polynomial<ProductElement<T>>, Integer> result = new TreeMap<>();
		for (Polynomial<T> prime : factorization.primeFactors()) {
			result.put(getUnivariatePolynomialRing().getEmbedding(prime, injection(nonZero.get())),
					factorization.multiplicity(prime));
		}
		return new FactorizationResult<>(inject(factorization.getUnit(), nonZero.get()), result);
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		for (R ring : rings) {
			if (!ring.isPrincipalIdealDomain()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isDedekindDomain() {
		Optional<R> nonZero = uniqueNonZeroRing();
		return nonZero.isPresent() && nonZero.get().isDedekindDomain();
	}

	@Override
	public boolean isDivisible(ProductElement<T> dividend, ProductElement<T> divisor) {
		for (int i = 0; i < numberOfFactors(); i++) {
			if (!rings.get(i).isDivisible(dividend.get(i), divisor.get(i))) {
				return false;
			}
		}
		return true;
	}

	@Override
	public QuotientAndRemainderResult<ProductElement<T>> quotientAndRemainder(ProductElement<T> dividend,
			ProductElement<T> divisor) {
		List<T> quotientValues = new ArrayList<>();
		List<T> remainderValues = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			QuotientAndRemainderResult<T> qr = rings.get(i).quotientAndRemainder(dividend.get(i), divisor.get(i));
			quotientValues.add(qr.getQuotient());
			remainderValues.add(qr.getRemainder());
		}
		return new QuotientAndRemainderResult<>(getElement(quotientValues), getElement(remainderValues));
	}

	@Override
	public BigInteger euclidMeasure(ProductElement<T> t) {
		Optional<Integer> nonZero = uniqueNonZeroRingIndex();
		if (nonZero.isEmpty()) {
			return null;
		}
		return rings.get(nonZero.get()).euclidMeasure(t.get(nonZero.get()));
	}

	@Override
	public ProductElement<T> projectToUnit(ProductElement<T> t) {
		List<T> units = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			units.add(rings.get(i).projectToUnit(t.get(i)));
		}
		return getElement(units);
	}

	@Override
	public Iterable<ProductElement<T>> getUnits() {
		List<Group<T>> units = new ArrayList<>();
		for (R ring : rings) {
			units.add(ring.getMultiplicativeGroup());
		}
		return new Iterable<ProductElement<T>>() {
			@Override
			public Iterator<ProductElement<T>> iterator() {
				return ProductRing.this.iterator(units);
			}
		};
	}

	@Override
	public int krullDimension() {
		int dimension = -1;
		for (R ring : rings) {
			if (ring.krullDimension() > dimension) {
				dimension = ring.krullDimension();
			}
		}
		return dimension;
	}

	public ProductIdeal<T, R> getProductIdeal(List<Ideal<T>> ideals) {
		return new ProductIdeal<>(ideals, this);
	}

	public ProductIdeal<T, R> injectFillWithUnitIdeals(Ideal<T> ideal, int index) {
		List<Ideal<T>> ideals = new ArrayList<>();
		for (int i = 0; i < index; i++) {
			ideals.add(rings.get(i).getUnitIdeal());
		}
		ideals.add(ideal);
		for (int i = index + 1; i < numberOfFactors(); i++) {
			ideals.add(rings.get(i).getUnitIdeal());
		}
		return getProductIdeal(ideals);
	}

	@SuppressWarnings("unchecked")
	@Override
	public List<Ideal<ProductElement<T>>> maximalPrimeIdealChain(Ideal<ProductElement<T>> start) {
		ProductIdeal<T, R> ideal = (ProductIdeal<T, R>) start;
		List<Ideal<T>> maximum = Collections.emptyList();
		int index = -1;
		for (int i = 0; i < numberOfFactors(); i++) {
			List<Ideal<T>> chainThisFactor = rings.get(i).maximalPrimeIdealChain(ideal.projection(i));
			if (chainThisFactor.size() > maximum.size()) {
				maximum = chainThisFactor;
				index = i;
			}
		}
		List<Ideal<ProductElement<T>>> result = new ArrayList<>();
		for (Ideal<T> primeIdeal : maximum) {
			result.add(injectFillWithUnitIdeals(primeIdeal, index));
		}
		return result;
	}

	@SuppressWarnings("unchecked")
	@Override
	public List<Ideal<ProductElement<T>>> maximalPrimeIdealChain(Ideal<ProductElement<T>> start,
			Ideal<ProductElement<T>> end) {
		ProductIdeal<T, R> startIdeal = (ProductIdeal<T, R>) start;
		ProductIdeal<T, R> endIdeal = (ProductIdeal<T, R>) end;
		List<Ideal<T>> maximum = Collections.emptyList();
		int index = -1;
		for (int i = 0; i < numberOfFactors(); i++) {
			List<Ideal<T>> chainThisFactor = rings.get(i).maximalPrimeIdealChain(startIdeal.projection(i),
					endIdeal.projection(index));
			if (chainThisFactor.size() > maximum.size()) {
				maximum = chainThisFactor;
				index = i;
			}
		}
		List<Ideal<ProductElement<T>>> result = new ArrayList<>();
		for (Ideal<T> primeIdeal : maximum) {
			result.add(injectFillWithUnitIdeals(primeIdeal, index));
		}
		return result;
	}

	public List<List<T>> getColumns(List<ProductElement<T>> t) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			List<T> column = new ArrayList<>();
			for (ProductElement<T> s : t) {
				column.add(s.get(i));
			}
			result.add(column);
		}
		return result;
	}

	@Override
	public IdealResult<ProductElement<T>, ProductIdeal<T, R>> getIdealWithTransforms(
			List<ProductElement<T>> generators) {
		List<Ideal<T>> idealList = new ArrayList<>();
		List<List<List<T>>> transforms = new ArrayList<>();
		List<Vector<ProductElement<T>>> syzygies = new ArrayList<>();
		List<List<T>> generatorList = getColumns(generators);
		int maximumGeneratorSize = 0;
		for (int i = 0; i < numberOfFactors(); i++) {
			IdealResult<T, ? extends Ideal<T>> ideal = rings.get(i).getIdealWithTransforms(generatorList.get(i));
			idealList.add(ideal.getIdeal());
			transforms.add(ideal.getGeneratorExpressions());
			if (maximumGeneratorSize < ideal.getGeneratorExpressions().size()) {
				maximumGeneratorSize = ideal.getGeneratorExpressions().size();
			}
			for (Vector<T> syzygy : ideal.getSyzygies()) {
				List<ProductElement<T>> productSyzygy = new ArrayList<>();
				for (int j = 0; j < generators.size(); j++) {
					productSyzygy.add(inject(syzygy.get(j + 1), i));
				}
				syzygies.add(new Vector<>(productSyzygy));
			}
		}
		List<List<ProductElement<T>>> productTransforms = new ArrayList<>();
		for (int i = 0; i < maximumGeneratorSize; i++) {
			List<ProductElement<T>> generatorTransform = new ArrayList<>();
			for (int j = 0; j < generators.size(); j++) {
				List<T> values = new ArrayList<>();
				for (int k = 0; k < numberOfFactors(); k++) {
					values.add(transforms.get(k).get(i).get(j));
				}
				generatorTransform.add(getElement(values));
			}
			productTransforms.add(generatorTransform);
		}
		return new IdealResult<>(productTransforms, generators, getProductIdeal(idealList), syzygies);
	}

	@Override
	public ProductIdeal<T, R> getIdeal(List<ProductElement<T>> generators) {
		List<Ideal<T>> idealList = new ArrayList<>();
		List<List<T>> generatorList = getColumns(generators);
		for (int i = 0; i < numberOfFactors(); i++) {
			idealList.add(rings.get(i).getIdeal(generatorList.get(i)));
		}
		return getProductIdeal(idealList);
	}

	@SuppressWarnings("unchecked")
	@Override
	public ProductIdeal<T, R> intersect(Ideal<ProductElement<T>> t1, Ideal<ProductElement<T>> t2) {
		ProductIdeal<T, R> ideal1 = (ProductIdeal<T, R>) t1;
		ProductIdeal<T, R> ideal2 = (ProductIdeal<T, R>) t2;
		List<Ideal<T>> idealList = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			idealList.add(rings.get(i).intersect(ideal1.projection(i), ideal2.projection(i)));
		}
		return getProductIdeal(idealList);
	}

	@SuppressWarnings("unchecked")
	@Override
	public ProductIdeal<T, R> radical(Ideal<ProductElement<T>> t) {
		ProductIdeal<T, R> ideal = (ProductIdeal<T, R>) t;
		List<Ideal<T>> idealList = new ArrayList<>();
		for (int i = 0; i < numberOfFactors(); i++) {
			idealList.add(rings.get(i).radical(ideal.projection(i)));
		}
		return getProductIdeal(idealList);
	}

	@SuppressWarnings("unchecked")
	@Override
	public ModuloMaximalIdealResult<ProductElement<T>, ?, ProductRing<T, R>, ProductIdeal<T, R>, ?> moduloMaximalIdeal(
			Ideal<ProductElement<T>> ideal) {
		ProductIdeal<T, R> maximalIdeal = (ProductIdeal<T, R>) ideal;
		Optional<Integer> nonUnitIndex = maximalIdeal.uniqueNonUnitIdealIndex();
		if (nonUnitIndex.isEmpty()) {
			throw new ArithmeticException("Not a maximal ideal!");
		}
		return moduloMaximalIdeal(maximalIdeal, nonUnitIndex.get(),
				rings.get(nonUnitIndex.get()).moduloMaximalIdeal(maximalIdeal.projection(nonUnitIndex.get())));
	}

	private <S extends Element<S>, Ri extends Ring<T>, I extends Ideal<T>, F extends Field<S>> ModuloMaximalIdealResult<ProductElement<T>, S, ProductRing<T, R>, ProductIdeal<T, R>, F> moduloMaximalIdeal(
			ProductIdeal<T, R> ideal, int index, ModuloMaximalIdealResult<T, S, Ri, I, F> moduloRing) {
		return new ModuloMaximalIdealResult<>(this, ideal, moduloRing.getField(), new MathMap<>() {
			@Override
			public S evaluate(ProductElement<T> t) {
				return moduloRing.reduce(project(t, index));
			}
		}, new MathMap<>() {
			@Override
			public ProductElement<T> evaluate(S t) {
				return inject(moduloRing.lift(t), index);
			}
		});
	}

	public static class ProductIdeal<T extends Element<T>, R extends Ring<T>> extends AbstractIdeal<ProductElement<T>>
			implements Ideal<ProductElement<T>> {
		private List<Ideal<T>> ideals;
		private ProductRing<T, R> productRing;
		private List<ProductElement<T>> generators;

		private ProductIdeal(List<Ideal<T>> ideals, ProductRing<T, R> productRing) {
			super(productRing);
			this.ideals = ideals;
			this.productRing = productRing;
		}

		@Override
		public boolean isFinite() {
			for (Ideal<T> ideal : ideals) {
				if (!ideal.isFinite()) {
					return false;
				}
			}
			return true;
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			BigInteger number = BigInteger.ONE;
			for (Ideal<T> ideal : ideals) {
				number = number.multiply(ideal.getNumberOfElements());
			}
			return number;
		}

		private Optional<Integer> uniqueNonUnitIdealIndex() {
			int result = -1;
			for (int i = 0; i < ideals.size(); i++) {
				if (!ideals.get(i).contains(productRing.getRing(i).one())) {
					if (result != -1) {
						return Optional.empty();
					}
					result = i;
				}
			}
			return result == -1 ? Optional.empty() : Optional.of(result);
		}

		private Optional<Ideal<T>> uniqueNonUnitIdeal() {
			Optional<Integer> index = uniqueNonUnitIdealIndex();
			if (index.isEmpty()) {
				return Optional.empty();
			}
			return Optional.of(ideals.get(index.get()));
		}

		@Override
		public boolean isPrimary() {
			Optional<Ideal<T>> nonUnit = uniqueNonUnitIdeal();
			if (nonUnit.isEmpty()) {
				return false;
			}
			return nonUnit.get().isPrimary();
		}

		@Override
		public boolean isPrime() {
			Optional<Ideal<T>> nonUnit = uniqueNonUnitIdeal();
			if (nonUnit.isEmpty()) {
				return false;
			}
			return nonUnit.get().isPrime();
		}

		@Override
		public boolean isMaximal() {
			Optional<Ideal<T>> nonUnit = uniqueNonUnitIdeal();
			if (nonUnit.isEmpty()) {
				return false;
			}
			return nonUnit.get().isMaximal();
		}

		public Ideal<T> projection(int index) {
			return ideals.get(index);
		}

		private List<T> projectList(List<ProductElement<T>> list, int index) {
			List<T> result = new ArrayList<>();
			for (ProductElement<T> e : list) {
				result.add(productRing.project(e, index));
			}
			return result;
		}

		private List<ProductElement<T>> injectList(List<T> list, int index) {
			List<ProductElement<T>> result = new ArrayList<>();
			for (T e : list) {
				result.add(productRing.inject(e, index));
			}
			return result;
		}

//		@Override
//		public List<Vector<ProductElement<T>>> nonTrivialCombinations(List<ProductElement<T>> s) {
//			List<Vector<ProductElement<T>>> result = new ArrayList<>();
//			for (int i = 0; i < productRing.numberOfFactors(); i++) {
//				List<Vector<T>> factorResult = ideals.get(i).nonTrivialCombinations(projectList(s, i));
//				for (Vector<T> row : factorResult) {
//					result.add(new Vector<>(injectList(row.asList(), i)));
//				}
//			}
//			return result;
//		}

		@Override
		public List<ProductElement<T>> generators() {
			if (generators == null) {
				generators = new ArrayList<>();
				List<Iterator<T>> generatorIterator = new ArrayList<>();
				for (int i = 0; i < productRing.numberOfFactors(); i++) {
					generatorIterator.add(ideals.get(i).generators().iterator());
				}
				while (true) {
					boolean noNext = true;
					List<T> values = new ArrayList<>();
					for (int i = 0; i < productRing.numberOfFactors(); i++) {
						if (generatorIterator.get(i).hasNext()) {
							noNext = false;
							values.add(generatorIterator.get(i).next());
						} else {
							values.add(productRing.getRing(i).zero());
						}
					}
					if (noNext) {
						break;
					}
					generators.add(productRing.getElement(values));
				}
			}
			return generators;
		}

		@Override
		public List<ProductElement<T>> generate(ProductElement<T> t) {
			List<ProductElement<T>> result = new ArrayList<>();
			List<Iterator<T>> generateIterator = new ArrayList<>();
			for (int i = 0; i < productRing.numberOfFactors(); i++) {
				generateIterator.add(ideals.get(i).generate(t.get(i)).iterator());
			}
			while (true) {
				boolean noNext = true;
				List<T> values = new ArrayList<>();
				for (int i = 0; i < productRing.numberOfFactors(); i++) {
					if (generateIterator.get(i).hasNext()) {
						noNext = false;
						values.add(generateIterator.get(i).next());
					} else {
						values.add(productRing.getRing(i).zero());
					}
				}
				if (noNext) {
					break;
				}
				result.add(productRing.getElement(values));
			}
			return result;
		}

		@Override
		public ProductElement<T> residue(ProductElement<T> t) {
			List<T> values = new ArrayList<>();
			for (int i = 0; i < productRing.numberOfFactors(); i++) {
				values.add(ideals.get(i).residue(t.get(i)));
			}
			return productRing.getElement(values);
		}

		@Override
		public boolean contains(ProductElement<T> t) {
			for (int i = 0; i < productRing.numberOfFactors(); i++) {
				if (!ideals.get(i).contains(t.get(i))) {
					return false;
				}
			}
			return true;
		}

	}
}
