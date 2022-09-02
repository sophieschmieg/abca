package fields.finitefields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.finitefields.ModuloIntegerRing.ModuloIntegerRingElement;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractIdeal;
import fields.helper.AbstractRing;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.IntegerIdeal;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.local.Value;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class ModuloIntegerRing extends AbstractRing<ModuloIntegerRingElement> {
	private BigInteger n;

	public ModuloIntegerRing(int n) {
		this.n = BigInteger.valueOf(n);
	}

	public ModuloIntegerRing(BigInteger n) {
		this.n = n;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public ModuloIntegerRingElement zero() {
		return this.getElement(0);
	}

	@Override
	public ModuloIntegerRingElement one() {
		return this.getElement(1);
	}

	@Override
	public BigInteger characteristic() {
		return n;
	}

	@Override
	public ModuloIntegerRingElement add(ModuloIntegerRingElement t1, ModuloIntegerRingElement t2) {
		return getElement(t1.value.add(t2.value));
	}

	@Override
	public ModuloIntegerRingElement negative(ModuloIntegerRingElement t) {
		return this.getElement(n.subtract(t.value));
	}

	@Override
	public ModuloIntegerRingElement multiply(ModuloIntegerRingElement t1, ModuloIntegerRingElement t2) {
		return getElement(t1.value.multiply(t2.value));
	}

	@Override
	public ModuloIntegerRingElement inverse(ModuloIntegerRingElement t) {
		return getElement(t.value.modInverse(n));
	}

	@Override
	public ModuloIntegerRingElement getRandomElement() {
		return getElement(MiscAlgorithms.randomBigInteger(new Random(), n));
	}

	@Override
	public BigInteger getNumberOfElements() {
		return this.n;
	}

	public ModuloIntegerRingElement reduce(IntE t) {
		return new ModuloIntegerRingElement(t.getValue());
	}

	public MathMap<IntE, ModuloIntegerRingElement> reduce() {
		return new MathMap<>() {
			@Override
			public ModuloIntegerRingElement evaluate(IntE t) {
				return reduce(t);
			}
		};
	}

	public IntE lift(ModuloIntegerRingElement t) {
		return Integers.z().getInteger(t.getValue());
	}

	public MathMap<ModuloIntegerRingElement, IntE> lift() {
		return new MathMap<>() {
			@Override
			public IntE evaluate(ModuloIntegerRingElement t) {
				return lift(t);
			}
		};
	}

	public Optional<Fraction> rationalReconstruction(ModuloIntegerRingElement t, BigInteger numeratorBound,
			BigInteger denominatorBound) {
		return Rationals.q().rationalReconstruction(lift(t), new IntE(n), numeratorBound, denominatorBound);
	}

	public Optional<Fraction> rationalReconstruction(ModuloIntegerRingElement t) {
		return Rationals.q().rationalReconstruction(lift(t), new IntE(n));
	}

	public IntE getOrder(ModuloIntegerRingElement t) {
		Integers z = Integers.z();
		IntE order = z.eulerToitent(z.getInteger(n));
		FactorizationResult<IntE, IntE> toitentFactors = z.uniqueFactorization(order);
		for (IntE prime : toitentFactors.primeFactors()) {
			int power = toitentFactors.multiplicity(prime);
			while (power > 0 && power(t, z.divideChecked(order, prime)).equals(one())) {
				power--;
				order = z.divideChecked(order, prime);
			}
		}
		return order;
	}

	@Override
	public Iterator<ModuloIntegerRingElement> iterator() {
		return new Iterator<ModuloIntegerRingElement>() {
			private BigInteger i = BigInteger.ZERO;

			@Override
			public boolean hasNext() {
				return this.i.compareTo(n) < 0;
			}

			@Override
			public ModuloIntegerRingElement next() {
				this.i = this.i.add(BigInteger.ONE);
				return ModuloIntegerRing.this.getElement(i.subtract(BigInteger.ONE));
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public String toString() {
		return "Z/" + n + "Z";
	}

	public ModuloIntegerRingElement getElement(int value) {
		return new ModuloIntegerRingElement(BigInteger.valueOf(value));
	}

	public ModuloIntegerRingElement getElement(BigInteger value) {
		return new ModuloIntegerRingElement(value);
	}

	@Override
	public boolean isUnit(ModuloIntegerRingElement t) {
		return n.gcd(t.value).equals(BigInteger.ONE);
	}

	@Override
	public BigInteger getNumberOfUnits() {
		Integers z = Integers.z();
		return z.eulerToitent(new IntE(n)).getValue();
	}

	@Override
	public Iterable<ModuloIntegerRingElement> getUnits() {
		return new Iterable<ModuloIntegerRingElement>() {

			@Override
			public Iterator<ModuloIntegerRingElement> iterator() {
				return new Iterator<ModuloIntegerRingElement>() {
					private BigInteger i = BigInteger.ZERO;

					@Override
					public boolean hasNext() {
						return this.i.compareTo(n) < 0;
					}

					@Override
					public ModuloIntegerRingElement next() {
						do {
							this.i = this.i.add(BigInteger.ONE);
						} while (!ModuloIntegerRing.this.isUnit(ModuloIntegerRing.this.getElement(i)));
						return ModuloIntegerRing.this.getElement(i);
					}
				};
			}
		};
	}

	@Override
	public boolean isIntegral() {
		return n.isProbablePrime(10);
	}

	@Override
	public boolean isReduced() {
		return Integers.z().uniqueFactorization(new IntE(n)).squareFree();
	}

	@Override
	public boolean isIrreducible() {
		return Integers.z().uniqueFactorization(new IntE(n)).primeFactors().size() == 1;
	}

	@Override
	public boolean isZeroDivisor(ModuloIntegerRingElement t) {
		return !isUnit(t);
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(ModuloIntegerRingElement dividend, ModuloIntegerRingElement divisor) {
		BigInteger divisorGcd = divisor.value.gcd(n);
		BigInteger dividendGcd = dividend.value.gcd(n);
		return dividendGcd.mod(divisorGcd).equals(BigInteger.ZERO);
	}

	@Override
	public QuotientAndRemainderResult<ModuloIntegerRingElement> quotientAndRemainder(ModuloIntegerRingElement dividend,
			ModuloIntegerRingElement divisor) {
		BigInteger divisorGcd = divisor.value.gcd(n);
		BigInteger dividendGcd = dividend.value.gcd(n);
		BigInteger[] qr = dividendGcd.divideAndRemainder(divisorGcd);
		if (!qr[1].equals(BigInteger.ZERO)) {
			throw new ArithmeticException("not divisble");
		}
		BigInteger reducedDividend = dividend.value.divide(divisorGcd);
		BigInteger reducedDivisor = divisor.value.divide(divisorGcd);
		BigInteger reducedN = n.divide(divisorGcd);
		ModuloIntegerRingElement result = new ModuloIntegerRingElement(
				reducedDividend.multiply(reducedDivisor.modInverse(reducedN)));
		return new QuotientAndRemainderResult<>(result, zero());
	}

	@Override
	public BigInteger euclidMeasure(ModuloIntegerRingElement t) {
		return null;
	}

	@Override
	public FactorizationResult<ModuloIntegerRingElement, ModuloIntegerRingElement> uniqueFactorization(
			ModuloIntegerRingElement t) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(z.getInteger(n));
		IntE e = z.getInteger(t.value);
		SortedMap<ModuloIntegerRingElement, Integer> result = new TreeMap<>();
		for (IntE prime : factors.primeFactors()) {
			int power = factors.multiplicity(prime);
			IntE factor = z.power(prime, power);
			IntE gcd = z.gcd(factor, e);
			if (!z.isUnit(gcd)) {
				e = z.divideChecked(e, gcd);
				int factorPower = 0;
				while (!z.isUnit(gcd)) {
					factorPower++;
					gcd = z.divideChecked(gcd, prime);
				}
				result.put(getInteger(prime), factorPower);
			}
		}
		if (!isUnit(getInteger(e))) {
			throw new ArithmeticException("Algorithm wrong");
		}
		return new FactorizationResult<>(getInteger(e), result);
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public ModuloIntegerRingElement projectToUnit(ModuloIntegerRingElement t) {
		if (t.equals(zero())) {
			return one();
		}
		return getElement(t.value.divide(t.value.gcd(n)));
	}

	public IdealResult<ModuloIntegerRingElement, ModularIntegerIdeal> getIdealWithTransforms(
			List<ModuloIntegerRingElement> generators) {
		if (generators.size() == 0) {
			return new IdealResult<>(Collections.singletonList(Collections.emptyList()), generators,
					new ModularIntegerIdeal(0), Collections.emptyList());
		}
		Integers z = Integers.z();
		List<IntE> lifted = new ArrayList<>();
		for (ModuloIntegerRingElement g : generators) {
			lifted.add(z.lift(g));
		}
		lifted.add(z.getInteger(n));
		IdealResult<IntE, IntegerIdeal> integerIdeal = z.getIdealWithTransforms(lifted);
		List<List<ModuloIntegerRingElement>> expressions = new ArrayList<>();
		for (List<IntE> integerList : integerIdeal.getGeneratorExpressions()) {
			List<ModuloIntegerRingElement> expression = new ArrayList<>();
			for (int i = 0; i < generators.size(); i++) {
				expression.add(reduce(integerList.get(i)));
			}
			expressions.add(expression);
		}
		List<Vector<ModuloIntegerRingElement>> syzygies = new ArrayList<>();
		for (Vector<IntE> integerSyzygy : integerIdeal.getSyzygies()) {
			List<ModuloIntegerRingElement> syzygy = new ArrayList<>();
			boolean nonZero = false;
			for (int i = 0; i < generators.size(); i++) {
				ModuloIntegerRingElement reduced = reduce(integerSyzygy.get(i + 1));
				syzygy.add(reduced);
				if (!reduced.equals(zero())) {
					nonZero = true;
				}
			}
			if (nonZero) {
				syzygies.add(new Vector<>(syzygy));
			}

		}
		return new IdealResult<>(expressions, generators,
				new ModularIntegerIdeal(reduce(integerIdeal.getIdeal().generators().get(0))), syzygies);
	}

	@Override
	public ModularIntegerIdeal getIdeal(List<ModuloIntegerRingElement> generators) {
		return getIdealWithTransforms(generators).getIdeal();
	}

	public ModularIntegerIdeal intersect(Ideal<ModuloIntegerRingElement> t1, Ideal<ModuloIntegerRingElement> t2) {
		return getIdeal(Collections.singletonList(lcm(t1.generators().get(0), t2.generators().get(0))));
	}

	@Override
	public ModularIntegerIdeal radical(Ideal<ModuloIntegerRingElement> t) {
		Integers z = Integers.z();
		IntE m = z.lift(t.generators().get(0));
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(m);
		ModuloIntegerRingElement radical = one();
		for (IntE prime : factors.primeFactors()) {
			radical = multiply(radical, reduce(prime));
		}
		return getIdeal(Collections.singletonList(radical));
	}

	@Override
	public PrimaryDecompositionResult<ModuloIntegerRingElement, ModularIntegerIdeal> primaryDecomposition(
			Ideal<ModuloIntegerRingElement> t) {
		Integers z = Integers.z();
		IntE generator = lift(t.generators().get(0));
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(generator);
		List<ModularIntegerIdeal> result = new ArrayList<>();
		List<ModularIntegerIdeal> radicals = new ArrayList<>();
		for (IntE prime : factors.primeFactors()) {
			result.add(getIdeal(Collections.singletonList(reduce(z.power(prime, factors.multiplicity(prime))))));
			radicals.add(getIdeal(Collections.singletonList(reduce(prime))));
		}
		return new PrimaryDecompositionResult<>(result, radicals);
	}

	@Override
	public ModuloMaximalIdealResult<ModuloIntegerRingElement, PFE, ModuloIntegerRing, ModularIntegerIdeal, PrimeField> moduloMaximalIdeal(
			Ideal<ModuloIntegerRingElement> ideal) {
		if (!ideal.isMaximal()) {
			throw new ArithmeticException("Not a maximal ideal!");
		}
		PrimeField fp = PrimeField.getPrimeField(ideal.generators().get(0).getValue());
		return new ModuloMaximalIdealResult<>(this, (ModularIntegerIdeal) ideal, fp, new MathMap<>() {

			@Override
			public PFE evaluate(ModuloIntegerRingElement t) {
				return fp.getElement(t.getValue());
			}
		}, new MathMap<>() {

			@Override
			public ModuloIntegerRingElement evaluate(PFE t) {
				return getElement(t.getValue());
			}
		});
	}

	@Override
	public ModuloIdealResult<ModuloIntegerRingElement, ?> moduloIdeal(Ideal<ModuloIntegerRingElement> ideal) {
		if (ideal.isPrime()) {
			ModuloMaximalIdealResult<ModuloIntegerRingElement, PFE, ModuloIntegerRing, ModularIntegerIdeal, PrimeField> mod = moduloMaximalIdeal(
					ideal);
			return new ModuloIdealResult<>(this, ideal, mod.getField(), mod.getReduction(), mod.getLift());
		}
		ModuloIntegerRing result = new ModuloIntegerRing(ideal.generators().get(0).getValue());
		return new ModuloIdealResult<>(this, ideal, result, new MathMap<>() {

			@Override
			public ModuloIntegerRingElement evaluate(ModuloIntegerRingElement t) {
				return result.getElement(t.getValue());
			}
		}, new MathMap<>() {
			@Override
			public ModuloIntegerRingElement evaluate(ModuloIntegerRingElement t) {
				return getElement(t.getValue());
			}
		});
	}

	@Override
	public Ideal<ModuloIntegerRingElement> getNilRadical() {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factorization = z.uniqueFactorization(z.getInteger(n));
		ModuloIntegerRingElement generator = one();
		for (IntE prime : factorization.primeFactors()) {
			generator = multiply(generator, reduce(prime));
		}
		return getIdeal(Collections.singletonList(generator));
	}

	public int krullDimension() {
		if (n.equals(BigInteger.ONE)) {
			return -1;
		}
		return 0;
	}

	@Override
	public List<Ideal<ModuloIntegerRingElement>> maximalPrimeIdealChain(Ideal<ModuloIntegerRingElement> start) {
		if (n.equals(BigInteger.ONE)) {
			return Collections.emptyList();
		}
		if (start.contains(one())) {
			throw new ArithmeticException("Not a proper ideal!");
		}
		PrimaryDecompositionResult<ModuloIntegerRingElement, ModularIntegerIdeal> decomposition = primaryDecomposition(
				start);
		return Collections.singletonList(decomposition.getRadicals().get(0));
	}

	@Override
	public List<Ideal<ModuloIntegerRingElement>> maximalPrimeIdealChain(Ideal<ModuloIntegerRingElement> start,
			Ideal<ModuloIntegerRingElement> end) {
		if (n.equals(BigInteger.ONE)) {
			return Collections.emptyList();
		}
		if (start.contains(one()) || end.contains(one())) {
			throw new ArithmeticException("Not a proper ideal!");
		}
		if (!end.contains(start)) {
			throw new ArithmeticException("Not an ideal chain");
		}
		PrimaryDecompositionResult<ModuloIntegerRingElement, ModularIntegerIdeal> decomposition = primaryDecomposition(
				start);
		return Collections.singletonList(decomposition.getRadicals().get(0));
	}

	public class ModuloIntegerRingElement extends AbstractElement<ModuloIntegerRingElement> {
		private BigInteger value;

		private ModuloIntegerRingElement(BigInteger value) {
			this.value = value.mod(n);
			while (this.value.compareTo(BigInteger.ZERO) < 0)
				this.value = this.value.add(n);
		}

		public BigInteger getValue() {
			return value;
		}

		public String toString() {
			return this.value.toString();
		}

		@Override
		public int compareTo(ModuloIntegerRingElement o) {
			return this.value.compareTo(o.value);
		}

	}

	public class ModularIntegerIdeal extends AbstractIdeal<ModuloIntegerRingElement> {
		private ModuloIntegerRingElement m;

		private ModularIntegerIdeal(ModuloIntegerRingElement m) {
			super(ModuloIntegerRing.this);
			this.m = m;
		}

		private ModularIntegerIdeal(int m) {
			this(BigInteger.valueOf(m));
		}

		private ModularIntegerIdeal(BigInteger m) {
			this(new ModuloIntegerRingElement(m));
		}

		@Override
		public boolean isFinite() {
			return true;
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			if (m.equals(zero())) {
				return BigInteger.ONE;
			}
			return n.divide(m.getValue());
		}

		@Override
		public Ideal<ModuloIntegerRingElement> annihilator() {
			if (m.equals(zero())) {
				getUnitIdeal();
			}
			return getIdeal(Collections.singletonList(getInteger(n.divide(m.getValue()))));
		}

//		@Override
//		public List<List<ModuloIntegerRingElement>> nonTrivialCombinations(List<ModuloIntegerRingElement> s) {
//			Integers z = Integers.z();
//			List<IntE> integerList = new ArrayList<>();
//			for (ModuloIntegerRingElement e : s) {
//				integerList.add(lift(e));
//			}
//			integerList.add(z.getInteger(n));
//			List<List<IntE>> integerResult = z.getUnitIdeal().nonTrivialCombinations(integerList);
//			List<List<ModuloIntegerRingElement>> result = new ArrayList<>();
//			for (List<IntE> integerRow : integerResult) {
//				List<ModuloIntegerRingElement> row = new ArrayList<>();
//				for (IntE value : integerRow.subList(0, s.size())) {
//					row.add(reduce(value));
//				}
//				result.add(row);
//			}
//			return result;
//		}

//		@Override
//		public List<Vector<ModuloIntegerRingElement>> getSyzygies() {
//			return Collections.singletonList(Collections.singletonList(getInteger(n.divide(m.getValue()))));
//		}

		@Override
		public List<ModuloIntegerRingElement> generators() {
			return Collections.singletonList(m);
		}

		@Override
		public List<ModuloIntegerRingElement> generate(ModuloIntegerRingElement t) {
			return Collections.singletonList(getElement(t.value.divide(m.value)));
		}

		@Override
		public ModuloIntegerRingElement residue(ModuloIntegerRingElement t) {
			return new ModuloIntegerRingElement(t.value.mod(m.value));
		}

		@Override
		public boolean isPrimary() {
			return Integers.z().getIdeal(Collections.singletonList(lift(m))).isPrimary();
		}

		@Override
		public boolean isPrime() {
			return isMaximal();
		}

		@Override
		public boolean isMaximal() {
			if (isIntegral()) {
				return m.equals(zero());
			}
			return m.value.isProbablePrime(100);
		}

		@Override
		public boolean contains(ModuloIntegerRingElement t) {
			return t.value.mod(m.value).equals(BigInteger.ZERO);
		}

		@Override
		public Value maximumPowerContains(ModuloIntegerRingElement t) {
			if (t.equals(zero()) || m.equals(one())) {
				return Value.INFINITY;
			}
			if (m.equals(zero())) {
				return Value.ZERO;
			}
			FactorizationResult<Ideal<ModuloIntegerRingElement>, Ideal<ModuloIntegerRingElement>> factors = idealFactorization(
					this);
			Value value = Value.INFINITY;
			Integers z = Integers.z();
			for (Ideal<ModuloIntegerRingElement> factor : factors.primeFactors()) {
				Ideal<IntE> liftedFactor = z.getIdeal(Collections.singletonList(z.lift(factor.generators().get(0))));
				value = value.min(new Value(
						z.valuation(z.lift(t), liftedFactor).value() / z.valuation(z.lift(m), liftedFactor).value()));
			}
			return value;
		}
	}

	private Matrix<IntE> liftAndAddModulus(Matrix<ModuloIntegerRingElement> m) {
		List<Vector<IntE>> lifted = new ArrayList<>();
		for (int i = 0; i < m.columns(); i++) {
			Vector<ModuloIntegerRingElement> column = m.column(i + 1);
			lifted.add(Vector.mapVector(lift(), column));
		}
		FreeModule<IntE> free = new FreeModule<>(Integers.z(), m.rows());
		for (int i = 0; i < m.rows(); i++) {
			lifted.add(free.scalarMultiply(n, free.getUnitVector(i + 1)));
		}
		return Matrix.fromColumns(lifted);
	}

	private Vector<IntE> liftVector(Vector<ModuloIntegerRingElement> t) {
		return Vector.mapVector(lift(), t);
	}

	private Vector<ModuloIntegerRingElement> reduceAndCutVector(Vector<IntE> t, int columns) {
		List<ModuloIntegerRingElement> result = new ArrayList<>();
		for (int i = 0; i < columns; i++) {
			result.add(reduce(t.get(i + 1)));
		}
		return new Vector<>(result);
	}

	@Override
	public boolean isSubModuleMember(MatrixModule<ModuloIntegerRingElement> module, Matrix<ModuloIntegerRingElement> m,
			Vector<ModuloIntegerRingElement> b) {
		Integers z = Integers.z();
		Matrix<IntE> liftedWithModulus = liftAndAddModulus(m);
		return z.isSubModuleMember(liftedWithModulus.getModule(z), liftedWithModulus, liftVector(b));
	}

	@Override
	public Vector<ModuloIntegerRingElement> asSubModuleMember(MatrixModule<ModuloIntegerRingElement> module,
			Matrix<ModuloIntegerRingElement> m, Vector<ModuloIntegerRingElement> b) {
		Integers z = Integers.z();
		Matrix<IntE> liftedWithModulus = liftAndAddModulus(m);
		return reduceAndCutVector(z.asSubModuleMember(liftedWithModulus.getModule(z), liftedWithModulus, liftVector(b)),
				m.columns());
	}

	@Override
	public List<Vector<ModuloIntegerRingElement>> syzygyProblem(MatrixModule<ModuloIntegerRingElement> module,
			Matrix<ModuloIntegerRingElement> m) {
		Integers z = Integers.z();
		Matrix<IntE> liftedWithModulus = liftAndAddModulus(m);
		List<Vector<ModuloIntegerRingElement>> result = new ArrayList<>();
		for (Vector<IntE> syzygy : z.syzygyProblem(liftedWithModulus.getModule(z), liftedWithModulus)) {
			Vector<ModuloIntegerRingElement> reduced = reduceAndCutVector(syzygy, m.columns());
			if (!reduced.equals(module.domain().zero())) {
				result.add(reduced);
			}
		}
		return result;
	}

	@Override
	public List<Vector<ModuloIntegerRingElement>> simplifySubModuleGenerators(
			MatrixModule<ModuloIntegerRingElement> module, Matrix<ModuloIntegerRingElement> m) {
		Integers z = Integers.z();
		Matrix<IntE> liftedWithModulus = liftAndAddModulus(m);
		List<Vector<ModuloIntegerRingElement>> result = new ArrayList<>();
		for (Vector<IntE> generator : z.simplifySubModuleGenerators(liftedWithModulus.getModule(z),
				liftedWithModulus)) {
			Vector<ModuloIntegerRingElement> reduced = Vector.mapVector(reduce(), generator);
			if (!reduced.equals(module.codomain().zero())) {
				result.add(reduced);
			}
		}
		return result;
	}
}
