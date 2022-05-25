package fields.helper;

import fields.helper.GenericExtensionField.GenericExtensionFieldElement;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.UnivariatePolynomial;

public class GenericExtensionField<T extends Element<T>>
		extends AbstractFieldExtension<T, GenericExtensionFieldElement<T>, GenericExtensionField<T>> {
//	private UnivariatePolynomial<T> minimalpolynomial;
//	private Field<T> field;
//	private GenericPrimeExtensionField<T> primeField;
//	private int degree;
//	private UnivariatePolynomialRing<T> ring;
//	private boolean small;
//	private List<GenericExtensionFieldElement<T>> elements;
//	private Map<GenericExtensionFieldElement<T>, Integer> index;
//	private int[][] additionTable;
//	private int[][] multiplicationTable;
//	private int[] inversionTable;
//	private int[] negationTable;
//	private int zeroIndex;
//	private int oneIndex;
//	private int alphaIndex;
//	private GenericExtensionFieldElement<T> zero;
//	private GenericExtensionFieldElement<T> one;
//	private GenericExtensionFieldElement<T> alpha;
//	private Map<Integer, GenericExtensionFieldElement<T>> alphaPowers;
//	private boolean powersReady;
//
	public GenericExtensionField(UnivariatePolynomial<T> minimalpolynomial, Field<T> field) {
		super(minimalpolynomial, field, "x");
	}

	@Override
	protected GenericExtensionFieldElement<T> fromSmallDegreePolynomial(UnivariatePolynomial<T> polynomial) {
		return new GenericExtensionFieldElement<>(this.getBaseField(), polynomial);
	}

	@Override
	public GenericExtensionField<T> makeExtension(UnivariatePolynomial<T> minimalPolynomial) {
		return new GenericExtensionField<>(minimalPolynomial, this.getBaseField());
	}

//
//	public GenericExtensionField(UnivariatePolynomial<T> minimalpolynomial, Field<T> field, boolean checkIrreducible,
//			boolean buildCache) {
//		super(new GenericPrimeExtensionField<>(field), field, minimalpolynomial.degree(), 1);
//		this.field = field;
//		this.powersReady = false;
//		this.primeField = (GenericPrimeExtensionField<T>) getBaseField();
//		this.ring = field.getUnivariatePolynomialRing();
//		this.minimalpolynomial = ring.normalize(minimalpolynomial);
//		if (checkIrreducible && field.isIrreducible(this.minimalpolynomial)) {
//			throw new ArithmeticException("Not irreducible!");
//		}
//		this.degree = this.minimalpolynomial.degree();
//		this.small = false;
//		this.alphaPowers = new TreeMap<>();
//		this.zero = getPrimeEmbedding(field.zero());
//		this.one = getPrimeEmbedding(field.one());
//		if (this.minimalpolynomial.degree() > 1) {
//			this.alpha = fromPolynomial(ring.getVar());
//		} else {
//			this.alpha = fromPolynomial(
//					ring.getEmbedding(field.negative(this.minimalpolynomial.univariateCoefficient(0))));
//		}
//		GenericExtensionFieldElement<T> power = one;
//		for (int i = 0; i < 2 * degree(); i++) {
//			this.alphaPowers.put(i, power);
//			power = multiply(power, alpha);
//		}
//		this.powersReady = true;
//		if (buildCache && this.isFinite() && this.getNumberOfElements().compareTo(BigInteger.valueOf(512)) <= 0) {
//			int num = this.getNumberOfElements().intValueExact();
//			this.elements = new ArrayList<>();
//			this.index = new TreeMap<>();
//			this.additionTable = new int[num][num];
//			this.multiplicationTable = new int[num][num];
//			this.negationTable = new int[num];
//			this.inversionTable = new int[num];
//			FreeModule<T> free = new FreeModule<>(field, degree);
//			int i = 0;
//			for (Vector<T> element : free) {
//				UnivariatePolynomial<T> asPolynomial = ring.getPolynomial(element.asList());
//				GenericExtensionFieldElement<T> e = new GenericExtensionFieldElement<>(field, primeField,
//						minimalpolynomial, asPolynomial, i);
//				this.elements.add(e);
//				this.index.put(e, i);
//				i++;
//			}
//			i = 0;
//			for (GenericExtensionFieldElement<T> e : this.elements) {
//				int j = 0;
//				for (GenericExtensionFieldElement<T> f : this.elements) {
//					this.additionTable[i][j] = this.index.get(this.add(e, f));
//					this.multiplicationTable[i][j] = this.index.get(this.multiply(e, f));
//					j++;
//				}
//				if (e.equals(zero())) {
//					this.zeroIndex = i;
//				} else if (e.equals(one())) {
//					this.oneIndex = i;
//				} else if (e.equals(alpha())) {
//					this.alphaIndex = i;
//				}
//				this.negationTable[i] = this.index.get(this.negative(e));
//				if (i != zeroIndex) {
//					this.inversionTable[i] = this.index.get(this.inverse(e));
//				}
//				i++;
//			}
//			this.small = true;
//		}
//	}
//
//	@Override
//	public Exactness exactness() {
//		return field.exactness();
//	}
//
//	@Override
//	public Extension<GenericExtensionFieldElement<T>, T, GenericExtensionFieldElement<T>> getExtension(
//			UnivariatePolynomial<GenericExtensionFieldElement<T>> minimalPolynomial) {
//		throw new UnsupportedOperationException("Too much recursion");
//	}
//
//	@Override
//	public String toString() {
//		if (this.minimalpolynomial.degree() == 1) {
//			return this.field.toString();
//		}
//		return this.field.toString() + "[X]/(" + this.minimalpolynomial.toString() + ")";
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> getPrimeEmbedding(T t) {
//		if (small) {
//			return fromPolynomial(ring.getEmbedding(t));
//		}
//		return new GenericExtensionFieldElement<T>(this.field, this.primeField, this.minimalpolynomial,
//				ring.getEmbedding(t));
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> getEmbedding(GenericExtensionFieldElement<T> t) {
//		return getPrimeEmbedding(t.asPrimeFieldElement());
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> zero() {
//		if (small) {
//			return elements.get(zeroIndex);
//		}
//		return zero;
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> one() {
//		if (small) {
//			return elements.get(oneIndex);
//		}
//		return one;
//	}
//
//	public GenericExtensionFieldElement<T> alpha() {
//		if (small) {
//			return elements.get(alphaIndex);
//		}
//		return alpha;
//	}
//
//	public GenericExtensionFieldElement<T> beta() {
//		return zero();
//	}
//
//	public GenericExtensionFieldElement<T> gamma() {
//		if (small) {
//			return elements.get(alphaIndex);
//		}
//		return fromPolynomial(ring.getVar(1));
//	}
//
//	public UnivariatePolynomial<GenericExtensionFieldElement<T>> minimalPolynomial() {
//		return primeField.getUnivariatePolynomialRing().getEmbedding(minimalpolynomial, new MathMap<>() {
//			@Override
//			public GenericExtensionFieldElement<T> evaluate(T t) {
//				return getPrimeEmbedding(t);
//			}
//		});
//	}
//
//	@Override
//	public UnivariatePolynomial<T> minimalPolynomialOverPrime() {
//		return minimalpolynomial;
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> add(GenericExtensionFieldElement<T> t1, GenericExtensionFieldElement<T> t2) {
//		if (small) {
//			return elements.get(additionTable[t1.index][t2.index]);
//		}
//		return fromPolynomial(ring.add(t1.asPolynomial, t2.asPolynomial));
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> negative(GenericExtensionFieldElement<T> t) {
//		if (small) {
//			return elements.get(negationTable[t.index]);
//		}
//		return fromPolynomial(ring.negative(t.asPolynomial));
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> multiply(GenericExtensionFieldElement<T> t1,
//			GenericExtensionFieldElement<T> t2) {
//		if (small) {
//			return elements.get(multiplicationTable[t1.index][t2.index]);
//		}
//		return fromPolynomial(this.ring.multiply(t1.asPolynomial, t2.asPolynomial));
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> inverse(GenericExtensionFieldElement<T> t) {
//		if (small) {
//			if (t.index == zeroIndex) {
//				throw new ArithmeticException("Division by zero");
//			}
//			return elements.get(inversionTable[t.index]);
//		}
//		UnivariatePolynomialRing.ExtendedResultantResult<T> eres = this.ring.extendedResultant(minimalpolynomial,
//				t.asPolynomial);
//		if (eres.getResultant().equals(field.zero())) {
//			throw new ArithmeticException("Division by zero: " + t.asPolynomial);
//		}
//		GenericExtensionFieldElement<T> inv = this
//				.fromPolynomial(this.ring.divideScalar(eres.getCoeff2(), eres.getGcd().univariateCoefficient(0)));
//		return inv;
//	}
//
//	public Polynomial<T> asPolynomial(GenericExtensionFieldElement<T> t) {
//		return t.asPolynomial;
//	}
//
//	@Override
//	public Vector<GenericExtensionFieldElement<T>> asVector(GenericExtensionFieldElement<T> t) {
//		List<T> asList = asPrimeVectorGamma(t).asList();
//		List<GenericExtensionFieldElement<T>> asVector = new ArrayList<>();
//		for (T e : asList) {
//			asVector.add(getPrimeEmbedding(e));
//		}
//		return new Vector<>(asVector);
//	}
//
//	@Override
//	public Vector<T> asPrimeVectorAlphaBetaGeneric(GenericExtensionFieldElement<GenericExtensionFieldElement<T>> t) {
//		return asPrimeVectorGamma(t.asPrimeFieldElement());
//	}
//
//	@Override
//	public Vector<T> asPrimeVectorAlphaBeta(GenericExtensionFieldElement<T> t) {
//		return asPrimeVectorGamma(t);
//	}
//
//	@Override
//	public Vector<T> asPrimeVectorGamma(GenericExtensionFieldElement<T> t) {
//		return ring.asVector(t.asPolynomial, degree - 1);
//	}
//
//	public GenericExtensionFieldElement<T> fromVectorOverPrime(Vector<T> t) {
//		return fromPolynomial(ring.getPolynomial(t.asList()));
//	}
//
//	@Override
//	public Matrix<T> alphaBetaToGammaBaseChange() {
//		return new FiniteVectorSpace<>(getPrimeField(), degree).matrixAlgebra().one();
//	}
//
//	@Override
//	public Matrix<T> gammaToAlphaBetaBaseChange() {
//		return new FiniteVectorSpace<>(getPrimeField(), degree).matrixAlgebra().one();
//	}
//
//	public List<GenericExtensionFieldElement<T>> conjugates(GenericExtensionFieldElement<T> t) {
//		throw new UnsupportedOperationException("Not implemented");
//	}
//
//	public GenericExtensionFieldElement<T> fromPolynomial(Polynomial<T> t) {
//		if (small) {
//			return this.elements.get(this.index.get(new GenericExtensionFieldElement<T>(this.field, this.primeField,
//					this.minimalpolynomial, ring.toUnivariate(t))));
//		}
//		if (powersReady && t.degree() >= degree && t.degree() < 2 * degree()) {
//			UnivariatePolynomial<T> p = ring.toUnivariate(t);
//			List<T> c = new ArrayList<>();
//			for (int i = 0; i < degree; i++) {
//				c.add(p.univariateCoefficient(i));
//			}
//			for (int i = degree; i <= t.degree(); i++) {
//				T coeff = p.univariateCoefficient(i);
//				UnivariatePolynomial<T> powerPolynomial = alphaPowers.get(i).asPolynomial;
//				for (int j = 0; j < degree; j++) {
//					c.set(j, field.add(c.get(j), field.multiply(coeff, powerPolynomial.univariateCoefficient(j))));
//				}
//			}
//			t = ring.getPolynomial(c);
//		}
//		return new GenericExtensionFieldElement<T>(this.field, this.primeField, this.minimalpolynomial,
//				ring.toUnivariate(t));
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> getRandomElement() {
//		return fromPolynomial(ring.getRandomElement(degree - 1));
//	}
//
//	@Override
//	public Iterator<GenericExtensionFieldElement<T>> iterator() {
//		return new Iterator<GenericExtensionFieldElement<T>>() {
//			private List<Iterator<T>> it = null;
//			private List<T> list = new ArrayList<T>();
//
//			private void init() {
//				if (this.it != null)
//					return;
//				this.it = new ArrayList<Iterator<T>>();
//				for (int i = 0; i < degree; i++) {
//					this.it.add(field.iterator());
//					if (i == 0)
//						this.list.add(null);
//					else
//						this.list.add(this.it.get(i).next());
//				}
//			}
//
//			@Override
//			public boolean hasNext() {
//				init();
//				for (int i = 0; i < degree; i++) {
//					if (this.it.get(i).hasNext())
//						return true;
//				}
//				return false;
//			}
//
//			@Override
//			public GenericExtensionFieldElement<T> next() {
//				init();
//				boolean broken = false;
//				for (int i = 0; i < degree; i++) {
//					if (this.it.get(i).hasNext()) {
//						this.list.set(i, this.it.get(i).next());
//						for (int j = 0; j < i; j++) {
//							this.it.set(j, field.iterator());
//							this.list.set(j, this.it.get(j).next());
//						}
//						broken = true;
//						break;
//					}
//				}
//				if (!broken)
//					throw new RuntimeException();
//				return fromPolynomial(ring.getPolynomial(list));
//			}
//		};
//	}
//
//	@Override
//	public GaloisGroup<T, GenericExtensionFieldElement<T>> galoisGroup() {
//		throw new UnsupportedOperationException("Not implemented!");
//	}
//
//	@Override
//	public GenericExtensionField<GenericExtensionFieldElement<T>> asGenericExtensionField() {
//		throw new UnsupportedOperationException("Too much recursion");
//	}
//
//	@Override
//	public GenericExtensionField<T> asGenericExtensionFieldOverPrime() {
//		return this;
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> fromGenericExtensionFieldElement(
//			GenericExtensionFieldElement<GenericExtensionFieldElement<T>> t) {
//		return fromPolynomial(ring.getEmbedding(t.asPolynomial, primeField.asPrimeFieldElementMap()));
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> fromGenericPrimeExtensionFieldElement(GenericExtensionFieldElement<T> t) {
//		return t;
//	}
//
//	@Override
//	public GenericExtensionFieldElement<GenericExtensionFieldElement<T>> asGenericExtensionFieldElement(
//			GenericExtensionFieldElement<T> s) {
//		return s.asGenericExtensionFieldElement();
//	}
//
//	@Override
//	public GenericExtensionFieldElement<T> asGenericPrimeExtensionFieldElement(GenericExtensionFieldElement<T> t) {
//		return t;
//	}
//
//	public CoordinateRing<T> asCoordinateRing() {
//		return new CoordinateRing<>(ring, ring.getIdeal(Collections.singletonList(minimalpolynomial)));
//	}
//
//	@Override
//	public T asPrimeFieldElement(GenericExtensionFieldElement<T> s) {
//		throw new UnsupportedOperationException("Not a base field");
//	}
//
	public static class GenericExtensionFieldElement<T extends Element<T>>
			extends AbstractElement<GenericExtensionFieldElement<T>>
			implements AlgebraicExtensionElement<T, GenericExtensionFieldElement<T>> {
//		private GenericPrimeExtensionField<T> primeField;
		private UnivariatePolynomial<T> asPolynomial;
//		private UnivariatePolynomial<T> minimalPolynomial;
//		private int degree;
//		private int index;

//		private GenericExtensionFieldElement(Field<T> base, GenericPrimeExtensionField<T> primeField,
//				UnivariatePolynomial<T> minimalpolynomial, UnivariatePolynomial<T> asPolynomial) {
//			this.base = base;
//			this.primeField = primeField;
//			this.minimalPolynomial = minimalpolynomial;
//			UnivariatePolynomialRing<T> ring = base.getUnivariatePolynomialRing();
//			if (asPolynomial.degree() < minimalpolynomial.degree()) {
//				this.asPolynomial = asPolynomial;
//			} else {
//				this.asPolynomial = ring.toUnivariate(ring.remainder(asPolynomial, minimalpolynomial));
//			}
//			this.degree = minimalpolynomial.degree();
//			this.index = -1;
//		}
//
//		private GenericExtensionFieldElement(Field<T> base, GenericPrimeExtensionField<T> primeField,
//				UnivariatePolynomial<T> minimalpolynomial, UnivariatePolynomial<T> asPolynomial, int index) {
//			this(base, primeField, minimalpolynomial, asPolynomial);
//			this.index = index;
//		}
		private GenericExtensionFieldElement(Field<T> base, UnivariatePolynomial<T> asPolynomial) {
			this.asPolynomial = asPolynomial;
		}

		@Override
		public String toString() {
			String[] elements = new String[] { "α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ",
					"ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω" };
			return asPolynomial.toString(elements[0], true);
		}
//
//		@Override
//		public String toString(String[] generators, int index) {
//			String generator = generators[index];
//			boolean first = true;
//			boolean onlyzero = true;
//			StringBuffer buf = new StringBuffer();
//			for (int i = 0; i <= asPolynomial.degree(); i++) {
//				T c = asPolynomial.univariateCoefficient(i);
//				String asString;
//				if (c instanceof ExtensionFieldElement<?, ?>) {
//					asString = ((ExtensionFieldElement<?, ?>) c).toString(generators, index + 1);
//				} else {
//					asString = c.toString();
//				}
//				if (asString.contains(" ")) {
//					asString = "(" + asString + ")";
//				}
//				if (!c.equals(this.base.zero())) {
//					onlyzero = false;
//					if (first)
//						first = false;
//					else
//						buf.append(" + ");
//					if (!c.equals(this.base.one()) || i == 0)
//						buf.append(asString);
//					if (i == 1)
//						buf.append(generator);
//					if (i > 1)
//						buf.append(generator + "^" + i);
//				}
//			}
//			if (onlyzero)
//				return "0";
//			return buf.toString();
//		}

		@Override
		public int compareTo(GenericExtensionFieldElement<T> o) {
			return asPolynomial.compareTo(o.asPolynomial);
		}

		public UnivariatePolynomial<T> asPolynomial() {
			return asPolynomial;
		}

//		@Override
//		public GenericExtensionFieldElement<GenericExtensionFieldElement<T>> asGenericExtensionFieldElement() {
//			UnivariatePolynomialRing<GenericExtensionFieldElement<T>> r = primeField.getUnivariatePolynomialRing();
//			return new GenericExtensionFieldElement<>(primeField, new GenericPrimeExtensionField<>(primeField),
//					r.getEmbedding(minimalPolynomial, primeField.getPrimeEmbeddingMap()),
//					r.getEmbedding(asPolynomial, primeField.getPrimeEmbeddingMap()));
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> asPrimeExtensionFieldElement() {
//			return this;
//		}
//
//		@Override
//		public T asPrimeFieldElement() {
//			return asPolynomial.univariateCoefficient(0);
//		}
	}
//
//	private static class GenericPrimeExtensionField<T extends Element<T>>
//			extends SkeletonFieldExtension<T, GenericExtensionFieldElement<T>> {
//		private Field<T> prime;
//		private UnivariatePolynomialRing<T> ring;
//
//		private GenericPrimeExtensionField(Field<T> prime) {
//			super(prime);
//			this.prime = prime;
//			this.ring = prime.getUnivariatePolynomialRing();
//		}
//
//		@Override
//		public Exactness exactness() {
//			return prime.exactness();
//		}
//
//		@Override
//		public GaloisGroup<T, GenericExtensionFieldElement<T>> galoisGroup() {
//			return new GaloisGroup<>(this, new FieldAutomorphism<>(this, new int[] { 0 }));
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> fromGenericPrimeExtensionFieldElement(
//				GenericExtensionFieldElement<T> t) {
//			return t;
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> asGenericPrimeExtensionFieldElement(GenericExtensionFieldElement<T> s) {
//			return s;
//		}
//
//		@Override
//		public List<GenericExtensionFieldElement<T>> conjugates(GenericExtensionFieldElement<T> s) {
//			return Collections.singletonList(s);
//		}
//
//		@Override
//		public Extension<GenericExtensionFieldElement<T>, T, GenericExtensionFieldElement<T>> getExtension(
//				UnivariatePolynomial<GenericExtensionFieldElement<T>> minimalPolynomial) {
//			throw new UnsupportedOperationException("Too much recursion");
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> getEmbedding(GenericExtensionFieldElement<T> t) {
//			return t;
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> getPrimeEmbedding(T t) {
//			return new GenericExtensionFieldElement<>(prime, this, ring.getVar(), ring.getEmbedding(t));
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> alpha() {
//			return zero();
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> beta() {
//			return zero();
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> gamma() {
//			return zero();
//		}
//
//		@Override
//		public UnivariatePolynomial<GenericExtensionFieldElement<T>> minimalPolynomial() {
//			return getUnivariatePolynomialRing().getVar();
//		}
//
//		@Override
//		public UnivariatePolynomial<T> minimalPolynomialOverPrime() {
//			return ring.getVar();
//		}
//
//		@Override
//		public GenericExtensionField<GenericExtensionFieldElement<T>> asGenericExtensionField() {
//			throw new UnsupportedOperationException("Too much recursion");
//		}
//
//		@Override
//		public GenericExtensionField<T> asGenericExtensionFieldOverPrime() {
//			throw new UnsupportedOperationException("Too much recursion");
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> fromGenericExtensionFieldElement(
//				GenericExtensionFieldElement<GenericExtensionFieldElement<T>> t) {
//			return t.asPrimeFieldElement();
//		}
//
//		@Override
//		public GenericExtensionFieldElement<GenericExtensionFieldElement<T>> asGenericExtensionFieldElement(
//				GenericExtensionFieldElement<T> s) {
//			throw new UnsupportedOperationException("Too much recursion");
//		}
//
//		@Override
//		public T asPrimeFieldElement(GenericExtensionFieldElement<T> s) {
//			return s.asPrimeFieldElement();
//		}
//
//		@Override
//		public Vector<GenericExtensionFieldElement<T>> asVector(GenericExtensionFieldElement<T> s) {
//			throw new UnsupportedOperationException("Too much recursion");
//		}
//
//		@Override
//		public Vector<T> asPrimeVectorAlphaBeta(GenericExtensionFieldElement<T> s) {
//			return new Vector<>(Collections.singletonList(asPrimeFieldElement(s)));
//		}
//
//		@Override
//		public Vector<T> asPrimeVectorAlphaBetaGeneric(
//				GenericExtensionFieldElement<GenericExtensionFieldElement<T>> s) {
//			return new Vector<>(Collections.singletonList(asPrimeFieldElement(s.asPrimeFieldElement())));
//		}
//
//		@Override
//		public Vector<T> asPrimeVectorGamma(GenericExtensionFieldElement<T> s) {
//			return new Vector<>(Collections.singletonList(asPrimeFieldElement(s)));
//		}
//
//		@Override
//		public Matrix<T> alphaBetaToGammaBaseChange() {
//			return asPrimeVectorSpace().matrixAlgebra().one();
//		}
//
//		@Override
//		public Matrix<T> gammaToAlphaBetaBaseChange() {
//			return asPrimeVectorSpace().matrixAlgebra().one();
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> zero() {
//			return getPrimeEmbedding(prime.zero());
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> one() {
//			return getPrimeEmbedding(prime.one());
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> add(GenericExtensionFieldElement<T> t1,
//				GenericExtensionFieldElement<T> t2) {
//			return getPrimeEmbedding(prime.add(t1.asPrimeFieldElement(), t2.asPrimeFieldElement()));
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> negative(GenericExtensionFieldElement<T> t) {
//			return getPrimeEmbedding(prime.negative(t.asPrimeFieldElement()));
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> multiply(GenericExtensionFieldElement<T> t1,
//				GenericExtensionFieldElement<T> t2) {
//			return getPrimeEmbedding(prime.multiply(t1.asPrimeFieldElement(), t2.asPrimeFieldElement()));
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> inverse(GenericExtensionFieldElement<T> t) {
//			return getPrimeEmbedding(prime.inverse(t.asPrimeFieldElement()));
//		}
//
//		@Override
//		public GenericExtensionFieldElement<T> getRandomElement() {
//			return getPrimeEmbedding(prime.getRandomElement());
//		}
//
//		@Override
//		public Iterator<GenericExtensionFieldElement<T>> iterator() {
//			return new Iterator<>() {
//				Iterator<T> it = prime.iterator();
//
//				@Override
//				public boolean hasNext() {
//					return it.hasNext();
//				}
//
//				@Override
//				public GenericExtensionFieldElement<T> next() {
//					return getPrimeEmbedding(it.next());
//				}
//			};
//		}
//
//	}

	@Override
	protected GenericExtensionField<T> asExtensionType() {
		return this;
	}
}
