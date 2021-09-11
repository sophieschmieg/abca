package fields.quaternions;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.local.PAdicField;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class RationalQuaternions extends AbstractQuaternions<Fraction> {
	private Set<IntE> ramifiedPrimes;
	private IntE discriminant;
	private boolean definite;

	private RationalQuaternions(Fraction a, Fraction b) {
		super(Rationals.q(), a, b);
		Integers z = Integers.z();
		IntE aUpToSquares = modSquares(a);
		IntE bUpToSquares = modSquares(b);
		this.ramifiedPrimes = new TreeSet<>();
		Set<IntE> potentiallyRamified = new TreeSet<>();
		potentiallyRamified.addAll(z.uniqueFactorization(aUpToSquares).primeFactors());
		potentiallyRamified.addAll(z.uniqueFactorization(bUpToSquares).primeFactors());
		potentiallyRamified.add(z.getInteger(2));
		Reals r = Reals.r(64);
		this.definite = RealQuaternions.quaternions(r, r.getFraction(a), r.getFraction(b)).isIntegral();
		this.discriminant = z.one();
		for (IntE prime : potentiallyRamified) {
			int accuracy = (int) Math.max(
					(Math.log(aUpToSquares.getValue().abs().doubleValue()) / Math.log(prime.getValue().doubleValue())),
					(Math.log(bUpToSquares.getValue().abs().doubleValue()) / Math.log(prime.getValue().doubleValue())))
					+ 1;
			PAdicField field = new PAdicField(prime.getValue(), accuracy);
			if (PAdicQuaternions.quaternions(field, field.getInteger(aUpToSquares), field.getInteger(bUpToSquares)).isIntegral()) {
				ramifiedPrimes.add(prime);
				this.discriminant = z.multiply(discriminant, prime);
			}
		}
		if ((ramifiedPrimes.size() % 2 == 0 && definite) || (ramifiedPrimes.size() % 2 == 1 && !definite)) {
			throw new ArithmeticException("parity wrong!");
		}
	}

	public static RationalQuaternions quaternions(IntE discriminant) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> discFactors = z.uniqueFactorization(discriminant);
		if (!discFactors.getUnit().equals(z.one())) {
			throw new ArithmeticException("Discriminant negative!");
		}
		Set<IntE> primes = discFactors.primeFactors();
		for (IntE prime : primes) {
			if (discFactors.multiplicity(prime) > 1) {
				throw new ArithmeticException("Discriminant not square free!");
			}
		}
		return quaternions(primes);
	}

	public static RationalQuaternions quaternions(Set<IntE> places) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		IntE unit = places.size() % 2 == 0 ? z.one() : z.getInteger(-1);
		IntE discriminant = z.one();
		for (IntE place : places) {
			discriminant = z.multiply(discriminant, place);
		}
		IntE adjustedDiscriminant = z.multiply(unit, discriminant);
		boolean two = places.contains(z.getInteger(2));
		Iterator<IntE> it = z.primes();
		while (true) {
			IntE prime = it.next();
			IntE adjustedPrime = z.multiply(unit, prime);
			int mod8 = adjustedPrime.getValue().mod(BigInteger.valueOf(8)).intValueExact();
			if (two && mod8 != 5) {
				continue;
			}
			if (!two && mod8 != 1) {
				continue;
			}
			boolean okay = true;
			for (IntE place : places) {
				if (place.equals(z.getInteger(2))) {
					continue;
				}
				if (MiscAlgorithms.jacobiSymbol(adjustedPrime.getValue(), place.getValue()) != -1) {
					okay = false;
					break;
				}
			}
			if (okay) {
				return new RationalQuaternions(q.getInteger(adjustedPrime), q.getInteger(adjustedDiscriminant));
			}
		}
	}

	public static RationalQuaternions quaternions(Fraction a, Fraction b) {
		return new RationalQuaternions(a, b);
	}

	private static IntE modSquares(Fraction t) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factorsNum = z.uniqueFactorization(t.getNumerator());
		FactorizationResult<IntE, IntE> factorsDenom = z.uniqueFactorization(t.getDenominator());
		IntE result = z.multiply(factorsNum.getUnit(), factorsDenom.getUnit());
		Set<IntE> primes = new TreeSet<>();
		primes.addAll(factorsNum.primeFactors());
		primes.addAll(factorsDenom.primeFactors());
		for (IntE prime : primes) {
			if ((factorsNum.multiplicity(prime) + factorsDenom.multiplicity(prime)) % 2 == 1) {
				result = z.multiply(result, prime);
			}
		}
		return result;
	}

	@Override
	public int hilbertSymbol() {
		return 1;
	}

	@Override
	public Fraction discriminant() {
		return Rationals.q().getInteger(discriminant);
	}

	@Override
	public boolean isIntegral() {
		return ramifiedPrimes.size() > 0;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof RationalQuaternions)) {
			return false;
		}
		RationalQuaternions o = (RationalQuaternions) obj;
		if (this.definite != o.definite) {
			return false;
		}
		return this.ramifiedPrimes.equals(o.ramifiedPrimes);
	}
	
	protected  Quaternion<Fraction> normalize( Quaternion<Fraction> t) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		IntE denom = z.one();
		IntE content = z.zero();
		for (Fraction value : t.asVector().asList()) {
			denom = z.lcm(denom, value.getDenominator());
			content = z.gcd(content, value.getNumerator());
		}
		List<Fraction> asList = new ArrayList<>();
		for (Fraction value : t.asVector().asList()) {
			asList.add(q.divide(q.multiply(denom, value), q.getInteger(content)));
		}
		return fromVector(new Vector<>(asList));
	}

	
	@Override
	protected Quaternion<Fraction> splittingElement() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		IntE aUpToSquares = modSquares(a());
		Fraction aSqrt = q.sqrt(q.divide(a(), q.getInteger(aUpToSquares))).keySet().iterator().next();
		IntE bUpToSquares = modSquares(b());
		NumberField nf = new NumberField(q.getUnivariatePolynomialRing()
				.getPolynomial(q.negative(q.getInteger(aUpToSquares)), q.zero(), q.one()));
		NumberFieldIntegers maximalOrder = nf.maximalOrder();
		FactorizationResult<IntE, IntE> bFactors = z.uniqueFactorization(bUpToSquares);
		NumberFieldIdeal ideal = maximalOrder.getUnitIdeal();
		for (IntE prime : bFactors.primeFactors()) {
			ideal = maximalOrder.multiply(ideal, maximalOrder.idealsOver(z.getIdeal(prime)).get(0));
		}
		if (ideal.generators().size() == 1) {
			NFE result = ideal.generators().get(0);
			return normalize(getElement(q.zero(), q.zero(), q.multiply(result.asPolynomial().univariateCoefficient(0), aSqrt),
					result.asPolynomial().univariateCoefficient(1)));
		}
		SortedSet<IntE> primes = new TreeSet<>();
		Iterator<IntE> primeIt = z.primes();
		IntE minkowskiBound = nf.minkowskiBound().roundDown();
		List<PFE> rhsList = new ArrayList<>();
		PrimeField f2 = PrimeField.getPrimeField(2);
		rhsList.add(f2.one());
		while (primeIt.hasNext()) {
			IntE prime = primeIt.next();
			if (!z.gcd(prime, bUpToSquares).equals(z.one())) {
				continue;
			}
			if (prime.compareTo(z.multiply(4, minkowskiBound)) > 0) {
				break;
			}
			primes.add(prime);
			rhsList.add(f2.zero());
		}
		Vector<PFE> rhs = new Vector<>(rhsList);
		List<NFE> candidates = new ArrayList<>();
		List<Vector<PFE>> columns = new ArrayList<>();
		while (true) {
			List<NFE> fromVector = new ArrayList<>();
			for (int i = 0; i < ideal.generators().size(); i++) {
				List<IntE> nfeFromVector = new ArrayList<>();
				for (int j = 0; j < nf.degree(); j++) {
					nfeFromVector.add(z.getRandomElement(minkowskiBound));
				}
				fromVector.add(maximalOrder.fromVector(new Vector<>(nfeFromVector)));
			}
			NFE idealElement = ideal.fromVector(new Vector<>(fromVector));
			IntE norm = nf.norm(idealElement).asInteger();
			norm = z.divideChecked(norm, bUpToSquares);
			Optional<FactorizationResult<IntE, IntE>> factorizationIfSmooth = z.uniqueFactorizationIfSmooth(norm, primes);
			if (factorizationIfSmooth.isEmpty()) {
				continue;
			}
			List<PFE> sqr = new ArrayList<>();
			sqr.add(f2.one());
			for (IntE prime : primes) {
				sqr.add(f2.getInteger(factorizationIfSmooth.get().multiplicity(prime)));
			}
			columns.add(new Vector<>(sqr));
			candidates.add(idealElement);
			Matrix<PFE> m = Matrix.fromColumns(columns);
			MatrixModule<PFE> matrices = new MatrixModule<>(f2, primes.size() + 1, candidates.size());
			if (matrices.hasSolution(m, rhs)) {
				NFE solution = nf.one();
				Vector<PFE> solved = matrices.solve(m, rhs);
				for (int i = 0; i < solved.dimension(); i++) {
					if (solved.get(i + 1).equals(f2.one())) {
						solution = nf.multiply(candidates.get(i), solution);
					}
				}
				IntE solvedNorm = nf.norm(solution).asInteger();
				solvedNorm = z.divide(solvedNorm, bUpToSquares);
				return normalize(getElement(q.zero(), q.zero(),
						q.multiply(solution.asPolynomial().univariateCoefficient(0), aSqrt),
						solution.asPolynomial().univariateCoefficient(1)));
			}
		}
	}
}
