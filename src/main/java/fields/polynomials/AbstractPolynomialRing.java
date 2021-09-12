package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.FormalPowerSeries;
import fields.local.FormalPowerSeries.PowerSeries;
import fields.local.Value;
import fields.vectors.FreeModule;
import fields.vectors.Vector;
import varieties.affine.AffinePoint;

public abstract class AbstractPolynomialRing<T extends Element<T>> extends AbstractAlgebra<T, Polynomial<T>>
		implements PolynomialRing<T> {
//	private Cache<Pair<Polynomial<T>, Polynomial<T>>, QuotientAndRemainderResult<Polynomial<T>>> quotientAndRemainderCache;

	protected AbstractPolynomialRing(boolean generateUnivariatePolynomialRing) {
		super(generateUnivariatePolynomialRing);
//		quotientAndRemainderCache = new Cache<>(1000);
	}

	public static <T extends Element<T>> PolynomialRing<T> getPolynomialRing(Ring<T> ring, int numvars,
			Comparator<Monomial> comparator) {
		if (numvars == 1) {
			return ring.getUnivariatePolynomialRing();
		}
		return new MultivariatePolynomialRing<T>(ring, numvars, comparator);
	}

	@Override
	public final BigInteger characteristic() {
		return getRing().characteristic();
	}

	@Override
	public boolean isCommutative() {
		return getRing().isCommutative();
	}

	@Override
	public final boolean isFinite() {
		return false;
	}

	@Override
	public final BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public final boolean isIntegral() {
		return getRing().isIntegral();
	}

	@Override
	public final boolean isReduced() {
		return getRing().isReduced();
	}

	@Override
	public boolean isIrreducible() {
		return getRing().isIrreducible();
	}

	@Override
	public final boolean isFree() {
		return true;
	}

	@Override
	public Polynomial<T> getEmbedding(Polynomial<T> t) {
		if (t.getPolynomialRing() == this) {
			return t;
		}
		int[] map = new int[numberOfVariables()];
		for (int i = 0; i < numberOfVariables(); i++) {
			map[i] = i;
		}
		return getEmbedding(t, map);
	}

	@Override
	public final Polynomial<T> getEmbeddingWithElimination(Polynomial<T> t, int shift) {
		if (this.numberOfVariables() != t.numberOfVariables() + shift)
			throw new ArithmeticException("Number of variables don't match!");
		int[] map = new int[t.numberOfVariables()];
		for (int i = 0; i < t.numberOfVariables(); i++)
			map[i] = i + shift;
		return this.getEmbedding(t, map);
	}

	@Override
	public Polynomial<T> projectToUnit(Polynomial<T> t) {
		return getEmbedding(getRing().projectToUnit(t.leadingCoefficient()));
	}

	@Override
	public Polynomial<T> upToUnit(Polynomial<T> t) {
		return divideScalar(t, getRing().projectToUnit(t.leadingCoefficient()));
	}

	@Override
	public Polynomial<T> multiply(T t, Monomial m, Polynomial<T> p) {
		return multiply(getEmbedding(t, m.exponents()), p);
	}

	@Override
	public MultivariatePolynomial<T> homogenize(Polynomial<T> t) {
		Map<Monomial, T> coeff = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			coeff.put(m.homogenizeMonomial(t.degree()), t.coefficient(m));
		}
		return new MultivariatePolynomial<>(coeff,
				new MultivariatePolynomialRing<>(getRing(), numberOfVariables() + 1, getComparator()));
	}

	@Override
	public Polynomial<T> homogenize(Polynomial<T> t, int coord) {
		Map<Monomial, T> coeff = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			coeff.put(m.homogenizeMonomial(t.degree(), coord), t.coefficient(m));
		}
		return getPolynomial(coeff);
	}

	@Override
	public Polynomial<T> dehomogenize(Polynomial<T> t, int coord) {
		if (coord < 1 || coord > this.numberOfVariables() || !isHomogeneous(t))
			throw new ArithmeticException("Not possible");
		List<T> values = new ArrayList<T>();
		for (int i = 0; i < this.numberOfVariables(); i++)
			if (i == coord - 1)
				values.add(this.getRing().one());
			else
				values.add(null);
		return this.partiallyEvaluate(t, values);
	}

	@SafeVarargs
	public final Polynomial<T> getLinear(T... coeffs) {
		return this.getLinear(Arrays.asList(coeffs));
	}

	@SafeVarargs
	public final T evaluate(Polynomial<T> t, T... values) {
		return this.evaluate(t, Arrays.asList(values));
	}

	@Override
	public int krullDimension() {
		return getRing().krullDimension() + numberOfVariables();
	}

	@Override
	public QuotientAndRemainderResult<Polynomial<T>> quotientAndRemainder(Polynomial<T> dividend,
			Polynomial<T> divisor) {
//		Pair<Polynomial<T>, Polynomial<T>> lookupKey = new Pair<>(dividend, divisor);
//		Optional<QuotientAndRemainderResult<Polynomial<T>>> cached = quotientAndRemainderCache.lookup(lookupKey);
//		if (cached.isPresent()) {
//			return cached.get();
//		}
		GeneralQuotientAndRemainderResult<T> gqr = generalQuotientAndRemainder(dividend,
				Collections.singletonList(divisor));
		new QuotientAndRemainderResult<>(gqr.getQuotients().get(0), gqr.getRemainder());
		QuotientAndRemainderResult<Polynomial<T>> result = new QuotientAndRemainderResult<>(gqr.getQuotients().get(0),
				gqr.getRemainder());
//		quotientAndRemainderCache.insert(new Pair<>(dividend, divisor), result);
		return result;
	}

	public static class ReduceAndExpressResult<T extends Element<T>> {
		private Polynomial<T> reduced;
		private List<Polynomial<T>> express;

		public ReduceAndExpressResult(Polynomial<T> reduced, List<Polynomial<T>> express) {
			this.reduced = reduced;
			this.express = express;
		}

		public Polynomial<T> getReduced() {
			return reduced;
		}

		public List<Polynomial<T>> getExpression() {
			return express;
		}

	}

	public ReduceAndExpressResult<T> reduceAndExpress(Polynomial<T> polynomial,
			List<Polynomial<T>> polynomialExpression, List<Polynomial<T>> basis,
			List<List<Polynomial<T>>> basisExpressions) {
		if (basis.isEmpty()) {
			Polynomial<T> unit = inverse(projectToUnit(polynomial));
			if (unit.equals(one())) {
				return new ReduceAndExpressResult<>(polynomial, polynomialExpression);
			}
			List<Polynomial<T>> expression = new ArrayList<>();
			for (Polynomial<T> coefficient : polynomialExpression) {
				expression.add(multiply(unit, coefficient));
			}
			return new ReduceAndExpressResult<>(multiply(unit, polynomial), expression);
		}
		GeneralQuotientAndRemainderResult<T> gqr = generalQuotientAndRemainder(polynomial, basis);
		// f = a_1 b_1 + ... a_n_b_n + r' = a_1 (c_11 g_1 + c_12 g2 + ... c_1m g_m) +
		// ... + r'
		// f = d_1 g_1 + ... + d_m g_m
		// r' = ur
		// r = u^-1((d_1-a_1c11-a_2c21-...)g_1 + ...)
		Polynomial<T> unit = inverse(projectToUnit(gqr.getRemainder()));
		List<Polynomial<T>> expression = new ArrayList<>();
		for (int j = 0; j < polynomialExpression.size(); j++) {
			expression.add(multiply(unit, polynomialExpression.get(j)));
		}
		for (int i = 0; i < basis.size(); i++) {
			Polynomial<T> quotient = gqr.getQuotients().get(i);
			List<Polynomial<T>> basisExpression = basisExpressions.get(i);
			for (int j = 0; j < expression.size(); j++) {
				expression.set(j, subtract(expression.get(j), multiply(unit, quotient, basisExpression.get(j))));
			}
		}
		return new ReduceAndExpressResult<>(multiply(unit, gqr.getRemainder()), expression);
	}

	@Override
	public Polynomial<T> reduce(Polynomial<T> polynomial, List<Polynomial<T>> basis) {
		return generalQuotientAndRemainder(polynomial, basis).getRemainder();
	}

	@Override
	public T evaluate(Polynomial<T> t, Vector<T> ts) {
		return evaluate(t, ts.asList());
	}

	@Override
	public T evaluate(Polynomial<T> t, AffinePoint<T> ts) {
		return evaluate(t, ts.getCoords());
	}

	@Override
	public PolynomialIdeal<T> getIdeal(List<Polynomial<T>> generators) {
		return new PolynomialIdeal<>(this, generators);
	}

	@Override
	public IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransforms(List<Polynomial<T>> generators) {
		GroebnerBasis<T> basis = buchberger(generators);
		return new IdealResult<>(basis.getExpression(), generators, new PolynomialIdeal<T>(this, basis.getBasis()));
	}

	@Override
	public PolynomialIdeal<T> getIdeal(@SuppressWarnings("unchecked") Polynomial<T>... generators) {
		return new PolynomialIdeal<T>(this, Arrays.asList(generators));
	}

	@Override
	public IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransforms(
			@SuppressWarnings("unchecked") Polynomial<T>... generators) {
		return getIdealWithTransforms(Arrays.asList(generators));
	}

	@Override
	public PolynomialIdeal<T> add(Ideal<Polynomial<T>> t1, Ideal<Polynomial<T>> t2) {
		List<Polynomial<T>> combined = new ArrayList<>();
		combined.addAll(t1.generators());
		combined.addAll(t2.generators());
		return getIdeal(combined);
	}

	@Override
	public PolynomialIdeal<T> intersect(Ideal<Polynomial<T>> t1, Ideal<Polynomial<T>> t2) {
		List<Polynomial<T>> generators = new ArrayList<>();
		PolynomialRing<T> ring = this.addVariableWithElimination(1);
		Polynomial<T> t = ring.getVar(1);
		Polynomial<T> omt = ring.subtract(ring.one(), t);
		for (Polynomial<T> f : t1.generators()) {
			generators.add(ring.multiply(t, ring.getEmbeddingWithElimination(f, 1)));
		}
		for (Polynomial<T> g : t2.generators()) {
			generators.add(ring.multiply(omt, ring.getEmbeddingWithElimination(g, 1)));
		}
		Ideal<Polynomial<T>> intersectionIdeal = ring.getIdeal(generators);
		List<Polynomial<T>> intersectionGenerators = new ArrayList<Polynomial<T>>();
		for (Polynomial<T> b : intersectionIdeal.generators()) {
			if (b.leadingMonomial().exponents()[0] == 0) {
				intersectionGenerators.add(getEmbeddingWithElimination(b, -1));
			}
		}
		return getIdeal(intersectionGenerators);
	}

	@Override
	public Polynomial<T> radical(Polynomial<T> t) {
		if (t.equals(zero())) {
			return t;
		}
		FactorizationResult<Polynomial<T>, T> squareFreeFactors = squareFreeFactorization(t);
		Polynomial<T> result = one();
		for (Polynomial<T> factor : squareFreeFactors.primeFactors()) {
			result = multiply(result, factor);
		}
		return result;
	}

	@Override
	public PolynomialIdeal<T> radical(Ideal<Polynomial<T>> t) {
		List<Polynomial<T>> radicalGenerators = new ArrayList<>();
		for (Polynomial<T> generator : t.generators()) {
			radicalGenerators.add(radical(generator));
		}
		return (PolynomialIdeal<T>) getIdeal(radicalGenerators);
	}

	@Override
	public PolynomialIdeal<T> getZeroIdeal() {
		return getIdeal(Collections.emptyList());
	}

	@Override
	public PolynomialIdeal<T> getUnitIdeal() {
		return getIdeal(Collections.singletonList(one()));
	}

//	@Override
//	public GroebnerBasis<T> buchberger(List<Polynomial<T>> generators) {
//		List<Polynomial<T>> groebner = new ArrayList<Polynomial<T>>();
//		List<List<Polynomial<T>>> expressions = new ArrayList<>(); 
//		for (Polynomial<T> basispolynomial : generators) {
//			Polynomial<T> reduced = this.reduce(basispolynomial, groebner);
//			if (!reduced.equals(this.zero())) {
//				groebner.add(upToUnit(reduced));
//			}
//		}
//		while (true) {
//			List<Polynomial<T>> newGroebner = new ArrayList<>();
//			for (int i = 1; i < groebner.size(); i++) {
//				for (int j = 0; j < i; j++) {
//					Polynomial<T> pi = groebner.get(i);
//					Polynomial<T> pj = groebner.get(j);
//					Monomial li = pi.leadingMonomial();
//					Monomial lj = pj.leadingMonomial();
//					T ci = pi.leadingCoefficient();
//					T cj = pj.leadingCoefficient();
//					T coefficient = getRing().lcm(pi.leadingCoefficient(), pj.leadingCoefficient());
//					Monomial lcm = li.lcm(lj);
//					li = lcm.divide(li);
//					lj = lcm.divide(lj);
//					ci = getRing().divide(coefficient, ci);
//					cj = getRing().divide(coefficient, cj);
//					Polynomial<T> pli = getEmbedding(ci, li.exponents());
//					Polynomial<T> plj = getEmbedding(cj, lj.exponents());
//					Polynomial<T> newGenerator = subtract(multiply(pli, pi), this.multiply(plj, pj));
//					newGenerator = reduce(newGenerator, groebner);
//					newGenerator = reduce(newGenerator, newGroebner);
//					if (!newGenerator.equals(this.zero())) {
//						newGroebner.add(upToUnit(newGenerator));
//					}
//					if (!getRing().isEuclidean() || getRing().krullDimension() == 0) {
//						continue;
//					}
//					Monomial libyj = li.divide(lj);
//					if (libyj != null) {
//						ExtendedEuclideanResult<T> extendedEuclid = getRing().extendedEuclidean(ci, cj);
//						if (!extendedEuclid.getGcd().equals(ci)) {
//							newGroebner.add(subtract(multiply(extendedEuclid.getCoeff1(), pi),
//									multiply(extendedEuclid.getCoeff2(), libyj, pj)));
//						}
//					}
//					Monomial ljbyi = lj.divide(li);
//					if (ljbyi != null) {
//						ExtendedEuclideanResult<T> extendedEuclid = getRing().extendedEuclidean(cj, ci);
//						if (!extendedEuclid.getGcd().equals(cj)) {
//							newGroebner.add(subtract(multiply(extendedEuclid.getCoeff1(), pj),
//									multiply(extendedEuclid.getCoeff2(), ljbyi, pi)));
//						}
//					}
//				}
//			}
//			newGroebner = reduceBasis(newGroebner);
//			newGroebner.addAll(groebner);
//			if (groebner.equals(newGroebner)) {
//				break;
//			}
//			groebner = newGroebner;
//		}
//		return reduceBasis(groebner);
//	}

	@Override
	public GroebnerBasis<T> buchberger(List<Polynomial<T>> generators) {
		List<Polynomial<T>> groebner = new ArrayList<>();
		List<List<Polynomial<T>>> expressions = new ArrayList<>();
		for (int i = 0; i < generators.size(); i++) {
			List<Polynomial<T>> generatorExpression = new ArrayList<>();
			for (List<Polynomial<T>> expression : expressions) {
				expression.add(zero());
			}
			for (int j = 0; j < i; j++) {
				generatorExpression.add(zero());
			}
			generatorExpression.add(one());
			ReduceAndExpressResult<T> reduced = this.reduceAndExpress(generators.get(i), generatorExpression, groebner,
					expressions);
			if (!reduced.getReduced().equals(this.zero())) {
				groebner.add(reduced.getReduced());
				expressions.add(reduced.getExpression());
			}
		}
		while (true) {
			List<Polynomial<T>> newGroebner = new ArrayList<>();
			List<List<Polynomial<T>>> newExpressions = new ArrayList<>();
			for (int i = 1; i < groebner.size(); i++) {
				for (int j = 0; j < i; j++) {
					Polynomial<T> pi = groebner.get(i);
					Polynomial<T> pj = groebner.get(j);
					Monomial li = pi.leadingMonomial();
					Monomial lj = pj.leadingMonomial();
					T ci = pi.leadingCoefficient();
					T cj = pj.leadingCoefficient();
					T coefficient = getRing().lcm(pi.leadingCoefficient(), pj.leadingCoefficient());
					Monomial lcm = li.lcm(lj);
					li = lcm.divide(li);
					lj = lcm.divide(lj);
					ci = getRing().divide(coefficient, ci);
					cj = getRing().divide(coefficient, cj);
					Polynomial<T> pli = getEmbedding(ci, li.exponents());
					Polynomial<T> plj = getEmbedding(cj, lj.exponents());
					Polynomial<T> newGenerator = subtract(multiply(pli, pi), this.multiply(plj, pj));
					List<Polynomial<T>> newGeneratorExpression = new ArrayList<>();
					for (int k = 0; k < generators.size(); k++) {
						newGeneratorExpression.add(subtract(multiply(pli, expressions.get(i).get(k)),
								this.multiply(plj, expressions.get(j).get(k))));
					}
					ReduceAndExpressResult<T> re = reduceAndExpress(newGenerator, newGeneratorExpression, groebner,
							expressions);
					newGenerator = re.getReduced();
					newGeneratorExpression = re.getExpression();
					re = reduceAndExpress(newGenerator, newGeneratorExpression, newGroebner, newExpressions);
					newGenerator = re.getReduced();
					newGeneratorExpression = re.getExpression();
					if (!newGenerator.equals(this.zero())) {
						newGroebner.add(newGenerator);
						newExpressions.add(newGeneratorExpression);
					}
					if (!getRing().isEuclidean() || getRing().krullDimension() == 0) {
						continue;
					}
					Monomial libyj = li.divide(lj);
					if (libyj != null) {
						ExtendedEuclideanResult<T> extendedEuclid = getRing().extendedEuclidean(ci, cj);
						if (!extendedEuclid.getGcd().equals(ci)) {
							newGroebner.add(subtract(multiply(extendedEuclid.getCoeff1(), pi),
									multiply(extendedEuclid.getCoeff2(), libyj, pj)));
							newGeneratorExpression = new ArrayList<>();
							for (int k = 0; k < generators.size(); k++) {
								newGeneratorExpression.add(subtract(
										multiply(extendedEuclid.getCoeff1(), expressions.get(i).get(k)),
										this.multiply(extendedEuclid.getCoeff2(), libyj, expressions.get(j).get(k))));

							}
							newExpressions.add(newGeneratorExpression);
						}
					}
					Monomial ljbyi = lj.divide(li);
					if (ljbyi != null) {
						ExtendedEuclideanResult<T> extendedEuclid = getRing().extendedEuclidean(cj, ci);
						if (!extendedEuclid.getGcd().equals(cj)) {
							newGroebner.add(subtract(multiply(extendedEuclid.getCoeff1(), pj),
									multiply(extendedEuclid.getCoeff2(), ljbyi, pi)));
							newGeneratorExpression = new ArrayList<>();
							for (int k = 0; k < generators.size(); k++) {
								newGeneratorExpression.add(subtract(
										multiply(extendedEuclid.getCoeff1(), expressions.get(j).get(k)),
										this.multiply(extendedEuclid.getCoeff2(), ljbyi, expressions.get(i).get(k))));

							}
							newExpressions.add(newGeneratorExpression);
						}
					}
				}
			}
			newGroebner.addAll(groebner);
			newExpressions.addAll(expressions);
			GroebnerBasis<T> b = reduceBasis(new GroebnerBasis<>(newGroebner, newExpressions));
			newGroebner = b.getBasis();
			newExpressions = b.getExpression();
			if (groebner.equals(newGroebner)) {
				break;
			}
			groebner = newGroebner;
			expressions = newExpressions;
		}
		return new GroebnerBasis<>(groebner, expressions);
	}

	@Override
	public GroebnerBasis<T> reduceBasis(GroebnerBasis<T> basis) {
		List<Polynomial<T>> list = new ArrayList<>();
		List<List<Polynomial<T>>> expressions = new ArrayList<>();
		list.addAll(basis.getBasis());
		expressions.addAll(basis.getExpression());
		SortedMap<Polynomial<T>, List<Polynomial<T>>> reduced = new TreeMap<>(Collections.reverseOrder());
		for (int i = 0; i < basis.getBasis().size(); i++) {
			list.remove(i);
			expressions.remove(i);
			ReduceAndExpressResult<T> reduce = this.reduceAndExpress(basis.getBasis().get(i),
					basis.getExpression().get(i), list, expressions);
			list.add(i, basis.getBasis().get(i));
			expressions.add(i, basis.getExpression().get(i));
			if (!reduce.getReduced().equals(this.zero())) {
				reduced.put(reduce.getReduced(), reduce.getExpression());
			}
		}
		list.clear();
		expressions.clear();
		list.addAll(reduced.keySet());
		expressions.addAll(reduced.values());
		return new GroebnerBasis<>(list, expressions);
	}

	@Override
	public boolean isGeneratingAlgebra(List<Polynomial<T>> s) {
		Ideal<Polynomial<T>> ideal = getIdeal(s);
		List<Polynomial<T>> variables = new ArrayList<>();
		for (int i = 0; i < numberOfVariables(); i++) {
			variables.add(getVar(i + 1));
		}
		Ideal<Polynomial<T>> all = getIdeal(variables);
		return ideal.equalsIdeal(all);
	}

	@Override
	public boolean isGeneratingModule(List<Polynomial<T>> s) {
		return false;
	}

	@Override
	public List<Polynomial<T>> getModuleGenerators() {
		throw new InfinityException();
	}

	@Override
	public String toString() {
		if (numberOfVariables() == 1) {
			return getRing().toString() + "[X]";
		} else if (numberOfVariables() == 2) {
			return getRing().toString() + "[X,Y]";
		} else if (numberOfVariables() == 3) {
			return getRing().toString() + "[X,Y,Z]";
		}
		StringBuffer buf = new StringBuffer();
		buf.append(getRing().toString() + "[");
		boolean first = true;
		for (int i = 0; i < numberOfVariables(); i++) {
			if (first) {
				first = false;
			} else {
				buf.append(",");
			}
			buf.append("X_" + i);
		}
		buf.append("]");
		return buf.toString();
	}

	@Override
	public final boolean hasCharacteristicRoot(Polynomial<T> t, int power) {
		if (characteristic().equals(BigInteger.ZERO)) {
			return false;
		}
		if (t.degree() <= 0) {
			return getRing().hasCharacteristicRoot(t.coefficient(getMonomial(new int[numberOfVariables()])), power);
		}
		for (Monomial m : t.monomials()) {
			if (t.coefficient(m).equals(getRing().zero())) {
				continue;
			}
			for (int i = 0; i < numberOfVariables(); i++) {
				if (!BigInteger.valueOf(m.exponents()[i]).mod(characteristic().pow(power)).equals(BigInteger.ZERO)) {
					return false;
				}
				if (!getRing().hasCharacteristicRoot(t.coefficient(m), power)) {
					return false;
				}
			}
		}
		return true;
	}

	@Override
	public final Polynomial<T> characteristicRoot(Polynomial<T> t, int power) {
		if (characteristic().equals(BigInteger.ZERO)) {
			throw new ArithmeticException("no 0th root!");
		}
		if (t.degree() == 0) {
			return getEmbedding(
					getRing().characteristicRoot(t.coefficient(getMonomial(new int[numberOfVariables()])), power));
		}
		Map<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			if (t.coefficient(m).equals(getRing().zero())) {
				continue;
			}
			int[] exp = new int[numberOfVariables()];
			for (int i = 0; i < numberOfVariables(); i++) {
				if (m.exponents()[i] % characteristic().pow(power).intValueExact() != 0) {
					throw new ArithmeticException("No pth root!");
				}
				exp[i] = m.exponents()[i] / characteristic().pow(power).intValueExact();
			}
			c.put(getMonomial(exp), getRing().characteristicRoot(t.coefficient(m), power));
		}
		return getPolynomial(c);
	}

	@Override
	public final boolean isEuclidean() {
		return this.numberOfVariables() == 1 && getRing().isIntegral() && getRing().krullDimension() == 0;
	}

	@Override
	public final boolean isPrincipalIdealDomain() {
		return isEuclidean();
	}

	@Override
	public final boolean isUniqueFactorizationDomain() {
		return getRing().isUniqueFactorizationDomain();
	}

	@Override
	public final boolean isDedekindDomain() {
		return isUniqueFactorizationDomain();
	}

	public Polynomial<T> moduloVarPower(Polynomial<T> t, int var, int power) {
		return moduloVarPowers(t, Collections.singletonMap(var, power));
	}

	public Polynomial<T> moduloVarPowers(Polynomial<T> t, Map<Integer, Integer> power) {
		SortedMap<Monomial, T> result = new TreeMap<>();
		int[] map = new int[numberOfVariables()];
		for (int i = 0; i < numberOfVariables(); i++) {
			if (power.containsKey(i + 1)) {
				map[i] = power.get(i + 1);
			} else {
				map[i] = -1;
			}
		}
		monomialLoop: for (Monomial m : t.monomials()) {
			for (int i = 0; i < numberOfVariables(); i++) {
				if (map[i] != -1 && m.exponents()[i] >= map[i]) {
					continue monomialLoop;
				}
			}
			result.put(m, t.coefficient(m));
		}
		return getPolynomial(result);
	}

	private List<Polynomial<T>> computeDeltas(List<Polynomial<T>> factors, List<Polynomial<T>> cofactors, int variable,
			List<Integer> degrees, Polynomial<T> target) {
		if (variable == 1) {
			UnivariatePolynomialRing<T> univariatePolynomialRing = getRing().getUnivariatePolynomialRing();
			List<Polynomial<T>> univariateCofactors = new ArrayList<>();
			for (Polynomial<T> cofactor : cofactors) {
				univariateCofactors.add(univariatePolynomialRing.getEmbedding(cofactor));
			}
			ExtendedEuclideanListResult<Polynomial<T>> ee = univariatePolynomialRing
					.extendedEuclidean(univariateCofactors);
			List<Polynomial<T>> result = new ArrayList<>();
			Polynomial<T> multiplier = multiply(target, inverse(ee.getGcd()));
			for (Polynomial<T> coefficient : ee.getCoeffs()) {
				result.add(multiply(multiplier, getEmbedding(coefficient)));
			}
			return result;
		}
		List<Polynomial<T>> reducedFactors = new ArrayList<>();
		List<Polynomial<T>> reducedCofactors = new ArrayList<>();
		Polynomial<T> reducedTarget = moduloVarPower(target, variable, 1);
		for (int i = 0; i < factors.size(); i++) {
			reducedFactors.add(moduloVarPower(factors.get(i), variable, 1));
			reducedCofactors.add(moduloVarPower(cofactors.get(i), variable, 1));
		}
		List<Polynomial<T>> deltas = computeDeltas(reducedFactors, reducedCofactors, variable - 1, degrees,
				reducedTarget);
		Map<Integer, Integer> degreeMap = new TreeMap<>();
		for (int i = 1; i < numberOfVariables(); i++) {
			degreeMap.put(i + 1, degrees.get(i - 1));
		}
		Polynomial<T> error = moduloVarPowers(target, degreeMap);
		for (int i = 0; i < cofactors.size(); i++) {
			error = moduloVarPowers(subtract(error, multiply(deltas.get(i), cofactors.get(i))), degreeMap);
		}
		for (int j = 1; j < degrees.get(variable - 2); j++) {
			List<Polynomial<T>> errorDelta = computeDeltas(reducedFactors, reducedCofactors, variable - 1, degrees,
					error);
			for (int i = 0; i < deltas.size(); i++) {
				Polynomial<T> diff = multiply(errorDelta.get(i), getVarPower(variable, j));
				deltas.set(i, add(deltas.get(i), diff));
				error = moduloVarPowers(subtract(error, multiply(diff, cofactors.get(i))), degreeMap);
			}
		}
		return deltas;
	}

	private Optional<List<Polynomial<T>>> henselLiftBivariatePolynomial(Polynomial<T> f, int accuracy) {
		Ring<T> ring = getRing();
		if (!ring.isEuclidean()) {
			throw new ArithmeticException("Base ring not euclidean");
		}
		if (numberOfVariables() != 2) {
			throw new ArithmeticException("Not bivariate");
		}
		UnivariatePolynomialRing<T> univariatePolynomialRing = ring.getUnivariatePolynomialRing();
		FormalPowerSeries<T> powerSeries = new FormalPowerSeries<>((Field<T>) ring, accuracy);
		UnivariatePolynomialRing<PowerSeries<T>> powerSeriesPolynomialRing = powerSeries.getUnivariatePolynomialRing();
		UnivariatePolynomialRing<Polynomial<T>> polynomialRingRing = ring.getUnivariatePolynomialRing()
				.getUnivariatePolynomialRing();
		UnivariatePolynomial<PowerSeries<T>> powerSeriesPolynomial = powerSeriesPolynomialRing
				.getEmbedding(asUnivariatePolynomial(f, 1), new MathMap<>() {

					@Override
					public PowerSeries<T> evaluate(Polynomial<T> t) {
						return powerSeries.getEmbedding(ring.getUnivariatePolynomialRing().toUnivariate(t));
					}
				});
		if (powerSeries.valuation(powerSeriesPolynomial.leadingCoefficient()).compareTo(Value.ZERO) > 0) {
			return Optional.empty();
		}
		powerSeriesPolynomial = powerSeriesPolynomialRing.normalize(powerSeriesPolynomial);
		f = fromUnivariatePolynomial(polynomialRingRing.getEmbedding(powerSeriesPolynomial, new MathMap<>() {
			@Override
			public Polynomial<T> evaluate(PowerSeries<T> t) {
				return powerSeries.roundToPolynomial(t, powerSeries.getAccuracy());
			}
		}), 1);
		Polynomial<T> modY = moduloVarPower(f, 2, 1);
		T unit = modY.leadingCoefficient();
		f = divideScalar(f, unit);
		modY = divideScalar(modY, unit);
		FactorizationResult<Polynomial<T>, T> modYFactors = ring
				.factorization(univariatePolynomialRing.toUnivariate(univariatePolynomialRing.getEmbedding(modY)));
		List<Polynomial<T>> factorList = new ArrayList<>();
		List<UnivariatePolynomial<T>> initialFactorList = new ArrayList<>();
		List<Polynomial<T>> cofactorList = new ArrayList<>();
		UnivariatePolynomial<T> initialProduct = univariatePolynomialRing.one();
		for (Polynomial<T> factor : modYFactors.primeFactors()) {
			Polynomial<T> embedded = getEmbedding(factor);
			initialProduct = univariatePolynomialRing.multiply(initialProduct, factor);
			if (modYFactors.multiplicity(factor) != 1) {
				return Optional.empty();
			}
			factorList.add(embedded);
			initialFactorList.add(univariatePolynomialRing.toUnivariate(factor));
		}
		List<UnivariatePolynomial<T>> delta = new ArrayList<>();
		for (Polynomial<T> factor : initialFactorList) {
			cofactorList.add(univariatePolynomialRing
					.toUnivariate(univariatePolynomialRing.divideChecked(initialProduct, factor)));
		}
		ExtendedEuclideanListResult<Polynomial<T>> ee = univariatePolynomialRing.extendedEuclidean(cofactorList);
		Polynomial<T> invertedGcd = univariatePolynomialRing.inverse(ee.getGcd());
		for (Polynomial<T> coeff : ee.getCoeffs()) {
			delta.add(univariatePolynomialRing.multiply(invertedGcd, coeff));
		}
		Polynomial<T> product = getEmbedding(initialProduct);
		for (int k = 1; k < accuracy; k++) {
			Polynomial<T> error = moduloVarPower(subtract(f, product), 2, k + 1);
			UnivariatePolynomial<T> univariateError = univariatePolynomialRing
					.toUnivariate(univariatePolynomialRing.getEmbedding(divideChecked(error, getVarPower(2, k))));
			product = one();
			for (int i = 0; i < delta.size(); i++) {
				Polynomial<T> factor = factorList.get(i);
				Polynomial<T> adjustedDelta = univariatePolynomialRing.multiply(delta.get(i), univariateError);
				adjustedDelta = univariatePolynomialRing.remainder(adjustedDelta, initialFactorList.get(i));
				factorList.set(i, add(multiply(getEmbedding(adjustedDelta), getVarPower(2, k)), factor));
				product = moduloVarPower(multiply(product, factorList.get(i)), 2, k + 2);
			}
		}
		return Optional.of(factorList);
	}

	private FactorizationResult<Polynomial<T>, T> henselLiftPolynomial(Polynomial<T> f, int variable, int accuracy) {
		Ring<T> ring = getRing();
		if (!ring.isEuclidean()) {
			throw new ArithmeticException("Base ring not euclidean");
		}
		UnivariatePolynomialRing<T> univariatePolynomialRing = ring.getUnivariatePolynomialRing();
		if (variable == 1) {
			FactorizationResult<Polynomial<T>, T> factorization = ring
					.factorization(univariatePolynomialRing.toUnivariate(univariatePolynomialRing.getEmbedding(f)));
			SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
			for (Polynomial<T> factor : factorization.primeFactors()) {
				result.put(getEmbedding(factor), factorization.multiplicity(factor));
			}
			return new FactorizationResult<>(factorization.getUnit(), result);
		}
		Polynomial<T> modY = moduloVarPower(f, variable, 1);
		FactorizationResult<Polynomial<T>, T> modYFactors = henselLiftPolynomial(modY, variable - 1, accuracy);
		List<Polynomial<T>> factorList = new ArrayList<>();
		List<Polynomial<T>> cofactorList = new ArrayList<>();
		Polynomial<T> product = getEmbedding(modYFactors.getUnit());
		for (Polynomial<T> factor : modYFactors.primeFactors()) {
			product = multiply(product, factor);
			if (modYFactors.multiplicity(factor) != 1) {
				throw new ArithmeticException("Not a square free polynomial!");
			}
			factorList.add(factor);
		}
		for (Polynomial<T> factor : factorList) {
			cofactorList.add(divideChecked(product, factor));
		}
		List<Polynomial<T>> initialFactorList = new ArrayList<>();
		initialFactorList.addAll(factorList);
		List<Polynomial<T>> delta = computeDeltas(factorList, cofactorList, variable, null, product);
		for (int k = 1; k < accuracy; k++) {
			Polynomial<T> error = moduloVarPower(subtract(f, product), 2, k + 1);
			product = getEmbedding(modYFactors.getUnit());
			for (int i = 0; i < delta.size(); i++) {
				Polynomial<T> factor = factorList.get(i);
				Polynomial<T> adjustedDelta = remainder(multiply(delta.get(i), error), initialFactorList.get(i));
				factorList.set(i, add(adjustedDelta, factor));
				product = multiply(product, factorList.get(i));
			}
		}
		SortedMap<Polynomial<T>, Integer> factorization = new TreeMap<>();
		for (Polynomial<T> factor : factorList) {
			factorization.put(factor, 1);
		}
		return new FactorizationResult<>(modYFactors.getUnit(), factorization);
	}

	private Iterator<Vector<T>> evaluationIterator(int dimension) {
		Ring<T> r = getRing();
		if (r.characteristic().equals(BigInteger.ZERO)) {
			return new Iterator<>() {
				private Iterator<UnivariatePolynomial<IntE>> it = Integers.z().getUnivariatePolynomialRing()
						.polynomials(dimension - 1);
				private UnivariatePolynomialRing<IntE> polynomials = Integers.z().getUnivariatePolynomialRing();

				@Override
				public boolean hasNext() {
					return true;
				}

				@Override
				public Vector<T> next() {
					Vector<IntE> next = polynomials.asVector(it.next(), dimension - 1);
					List<T> result = new ArrayList<>();
					for (IntE value : next.asList()) {
						result.add(r.getInteger(value));
					}
					return new Vector<>(result);
				}
			};
		}
		if (r.isFinite()) {
			return new FreeModule<>(r, dimension).iterator();
		}
		return new Iterator<>() {
			private FreeModule<T> module = new FreeModule<>(r, dimension);

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public Vector<T> next() {
				return module.getRandomElement();
			}
		};
	}

	private static class SquareFreeFactorizationResult<T extends Element<T>> {
		private List<Polynomial<T>> factors;
		private T unit;
	}

	private SquareFreeFactorizationResult<T> factorizeBivariateSquareFreePolynomial(Polynomial<T> t) {
		if (numberOfVariables() != 2) {
			throw new ArithmeticException("Not bivariate!");
		}
		if (!(getRing() instanceof Field<?>)) {
			throw new UnsupportedOperationException("Not yet implemented");
		}
		Ring<T> ring = getRing();
		UnivariatePolynomialRing<T> polynomialRing = ring.getUnivariatePolynomialRing();
		Iterator<Vector<T>> it = evaluationIterator(1);
		UnivariatePolynomial<Polynomial<T>> asUnivariate = asUnivariatePolynomial(t, 2);
		Polynomial<T> content = ring.getUnivariatePolynomialRing().getUnivariatePolynomialRing().content(asUnivariate);
		FactorizationResult<Polynomial<T>, T> contentFactors = ring
				.factorization(ring.getUnivariatePolynomialRing().toUnivariate(content));
		SquareFreeFactorizationResult<T> result = new SquareFreeFactorizationResult<>();
		result.unit = contentFactors.getUnit();
		result.factors = new ArrayList<>();
		for (Polynomial<T> contentFactor : contentFactors.primeFactors()) {
			result.factors.add(getEmbedding(contentFactor, new int[] { 0 }));
			if (contentFactors.multiplicity(contentFactor) != 1) {
				throw new ArithmeticException("Content not square free!");
			}
		}
		t = fromUnivariatePolynomial(
				ring.getUnivariatePolynomialRing().getUnivariatePolynomialRing().contentFree(asUnivariate), 2);
		asUnivariate = asUnivariatePolynomial(t, 1);
		content = ring.getUnivariatePolynomialRing().getUnivariatePolynomialRing().content(asUnivariate);
		contentFactors = ring.factorization(ring.getUnivariatePolynomialRing().toUnivariate(content));
		result.unit = ring.multiply(contentFactors.getUnit(), result.unit);
		for (Polynomial<T> contentFactor : contentFactors.primeFactors()) {
			result.factors.add(getEmbedding(contentFactor, new int[] { 1 }));
			if (contentFactors.multiplicity(contentFactor) != 1) {
				throw new ArithmeticException("Content not square free!");
			}
		}
		t = fromUnivariatePolynomial(
				ring.getUnivariatePolynomialRing().getUnivariatePolynomialRing().contentFree(asUnivariate), 1);
		if (t.degree() == 0) {
			result.unit = ring.multiply(result.unit, t.leadingCoefficient());
			return result;
		}
		int degreeX = t.degree(1);
		int degreeY = t.degree(2);
		while (it.hasNext()) {
			T point = it.next().get(1);
			List<T> evaluationPoint = new ArrayList<>();
			evaluationPoint.add(null);
			evaluationPoint.add(point);
			Polynomial<T> modYMinusPoint = partiallyEvaluate(t, evaluationPoint);
			if (modYMinusPoint.degree() != degreeX) {
				continue;
			}
			if (!polynomialRing
					.squareFreeFactorization(polynomialRing.toUnivariate(polynomialRing.getEmbedding(modYMinusPoint)))
					.squareFree()) {
				continue;
			}
			List<Polynomial<T>> substitute = new ArrayList<>();
			substitute.add(getVar(1));
			substitute.add(add(getVar(2), getEmbedding(point)));
			List<Polynomial<T>> backSubstitute = new ArrayList<>();
			backSubstitute.add(getVar(1));
			backSubstitute.add(subtract(getVar(2), getEmbedding(point)));
			Polynomial<T> substituted = substitute(t, substitute);
			int accuracy = 2 * degreeY + 2;
			Optional<List<Polynomial<T>>> maybeLiftedFactors = henselLiftBivariatePolynomial(substituted, accuracy);
			if (maybeLiftedFactors.isEmpty()) {
				continue;
			}
			List<Polynomial<T>> liftedFactors = maybeLiftedFactors.get();

//			Field<T> field = (Field<T>) getRing();
//			FormalPowerSeries<T> powerSeries = new FormalPowerSeries<>(field, 2 * t.degree(1) + 4);
//		UnivariatePolynomialRing<PowerSeries<T>> powerSeriesPolynomialRing = powerSeries.getUnivariatePolynomialRing();
//		UnivariatePolynomial<Polynomial<T>> asUnivariatePolynomial = asUnivariatePolynomial(t, 2);
//		UnivariatePolynomial<PowerSeries<T>> asPowerSeriesPolynomial = powerSeriesPolynomialRing
//				.getEmbedding(asUnivariatePolynomial, new MathMap<>() {
//					@Override
//					public PowerSeries<T> evaluate(Polynomial<T> c) {
//						return powerSeries.getEmbedding(field.getUnivariatePolynomialRing().toUnivariate(c));
//					}
//				});
//		asPowerSeriesPolynomial = powerSeriesPolynomialRing
//				.toUnivariate(powerSeriesPolynomialRing.normalize(asPowerSeriesPolynomial));
//		FactorizationResult<Polynomial<PowerSeries<T>>> factorization = powerSeries
//				.factorization(asPowerSeriesPolynomial);
//		List<Polynomial<PowerSeries<T>>> liftedFactors = new ArrayList<>();
//		for (Polynomial<PowerSeries<T>> factor : factorization.primeFactors()) {
//			if (factor.degree() > 0) {
//				liftedFactors.add(factor);
//			}
//		}
			Map<Integer, Polynomial<T>> powerSeriesFactors = new TreeMap<>();
			for (int i = 0; i < liftedFactors.size(); i++) {
				powerSeriesFactors.put(i, liftedFactors.get(i));
			}
			List<CombinedFactors<T>> powerSeriesCombinedFactors = new ArrayList<>();
			powerSeriesCombinedFactors.add(new CombinedFactors<T>(Collections.emptySortedSet(), one()));
			List<Polynomial<T>> polynomialFactors = new ArrayList<>();
			while (powerSeriesFactors.size() > 0) {
				CheckCombinationsResult<T> ccresult = checkCombinations(substituted, powerSeriesFactors,
						powerSeriesCombinedFactors, accuracy);
				substituted = ccresult.cofactor;
				for (Polynomial<T> factor : ccresult.factors) {
					polynomialFactors.add(substitute(factor, backSubstitute));
				}
				for (int k : ccresult.usedFactors) {
					powerSeriesFactors.remove(k);
				}
				powerSeriesCombinedFactors = ccresult.combined;
			}
			if (substituted.degree() != 0) {
				System.err.println(t);
				System.err.println(substituted);
				throw new ArithmeticException("recombination failed");
			}
			result.factors.addAll(polynomialFactors);
			result.unit = ring.multiply(result.unit, substituted.leadingCoefficient());
			return result;
		}
		// ran out of evaluation points, moving to field extension.
		return factorizeOverExtension(t);

	}

	private SquareFreeFactorizationResult<T> factorizeSquareFreePolynomial(Polynomial<T> t) {
		if (numberOfVariables() == 2) {
			return factorizeBivariateSquareFreePolynomial(t);
		}
		AbstractPolynomialRing<T> bivariate = (AbstractPolynomialRing<T>) AbstractPolynomialRing
				.getPolynomialRing(getRing(), 2, getComparator());
		UnivariatePolynomialRing<T> univariate = getRing().getUnivariatePolynomialRing();
		Iterator<Vector<T>> it = evaluationIterator(numberOfVariables() - 2);
		int degreeX = t.degree(1);
		if (degreeX == 0) {
			throw new ArithmeticException("Expected non trivial first variable!");
		}
		int degreeY = t.degree(2);
		Polynomial<T> leading = asUnivariatePolynomial(t, 1).leadingCoefficient();
		Polynomial<T> leadingRadical = eliminateVariable().radical(leading);
		mainLoop: while (it.hasNext()) {
			List<T> pointList = it.next().asList();
			for (int i = numberOfVariables(); i > 2; i--) {
				List<T> evaluationPoint = new ArrayList<>();
				for (int j = 0; j < i - 1; j++) {
					evaluationPoint.add(null);
				}
				evaluationPoint.addAll(pointList.subList(i - 3, pointList.size()));
				Polynomial<T> evaluated = bivariate.getEmbedding(partiallyEvaluate(t, evaluationPoint));
				if (evaluated.degree(i) != t.degree(i)) {
					continue mainLoop;
				}
			}
			List<T> evaluationPoint = new ArrayList<>();
			evaluationPoint.add(null);
			evaluationPoint.add(null);
			evaluationPoint.addAll(pointList);
			Polynomial<T> evaluated = bivariate.getEmbedding(partiallyEvaluate(t, evaluationPoint));
			if (evaluated.degree(1) != degreeX || evaluated.degree(2) != degreeY) {
				continue;
			}
			if (!bivariate.squareFreeFactorization(evaluated).squareFree()) {
				continue;
			}
			SquareFreeFactorizationResult<T> bivariateFactors = bivariate
					.factorizeBivariateSquareFreePolynomial(evaluated);
			Polynomial<T> leadingEvaluated = bivariate.asUnivariatePolynomial(evaluated, 1).leadingCoefficient();
			Polynomial<T> leadingEvaluatedRadical = univariate.radical(leadingEvaluated);
			Polynomial<T> leadingRadicalEvaluated = univariate.getEmbedding(
					eliminateVariable().evaluate(leading, evaluationPoint.subList(1, evaluationPoint.size())));
			Polynomial<T> m;
			List<Polynomial<T>> p = new ArrayList<>();
			if (!leadingEvaluatedRadical.equals(leadingRadicalEvaluated)) {
				m = leading;
				for (int i = 0; i < bivariateFactors.factors.size(); i++) {
					p.add(one());
				}
			} else {
				List<Polynomial<T>> squareFreeFactors = new ArrayList<>();
				for (Polynomial<T> bivariateFactor : bivariateFactors.factors) {
					Polynomial<T> leadingOfFactor = bivariate.asUnivariatePolynomial(bivariateFactor, 1)
							.leadingCoefficient();
					FactorizationResult<Polynomial<T>, T> squareFreeOfLeadingFactors = univariate
							.squareFreeFactorization(leadingOfFactor);
					squareFreeFactors.addAll(squareFreeOfLeadingFactors.primeFactors());
				}
				List<Polynomial<T>> gcdFreeBasis = gcdFreeBasis(squareFreeFactors);
			}
		}
		// ran out of evaluation points, moving to field extension.
		return factorizeOverExtension(t);
	}

	private List<Polynomial<T>> gcdFreeBasis(List<Polynomial<T>> set) {
		List<Polynomial<T>> reduced = new ArrayList<>();
		for (Polynomial<T> element : set) {
			if (element.degree() > 0) {
				reduced.add(element);
			}
		}
		int reducedSize = reduced.size();
		for (int i = 0; i < reducedSize; i++) {
			for (int j = 0; j < reducedSize; j++) {
				if (i == j) {
					continue;
				}
				Polynomial<T> gcd = gcd(reduced.get(i), reduced.get(j));
				if (gcd.degree() != 0) {
					reduced.add(gcd);
					reduced.set(i, divideChecked(reduced.get(i), gcd));
					reduced.set(j, divideChecked(reduced.get(j), gcd));
				}
			}
		}
		if (reducedSize == reduced.size()) {
			return reduced;
		}
		return gcdFreeBasis(reduced);
	}

	private SquareFreeFactorizationResult<T> factorizeOverExtension(Polynomial<T> t) {
		if (!getRing().isFinite() || !(getRing() instanceof Field<?>)) {
			throw new ArithmeticException("Should not run out of evaluation points in this case!");
		}
		Field<T> field = (Field<T>) getRing();
		T primitiveRoot = field.primitiveRoot();
		List<T> extensionCoefficients = new ArrayList<>();
		extensionCoefficients.add(field.negative(primitiveRoot));
		extensionCoefficients.add(field.zero());
		extensionCoefficients.add(field.one());
		return factorizeOverExtension(t,
				field.getExtension(field.getUnivariatePolynomialRing().getPolynomial(extensionCoefficients)));

	}

	private <Base extends Element<Base>, Ext extends AlgebraicExtensionElement<Base, Ext>, ExtField extends FieldExtension<Base, Ext, ExtField>> SquareFreeFactorizationResult<T> factorizeOverExtension(
			Polynomial<T> t, Extension<T, Base, Ext, ExtField> extension) {
		AbstractPolynomialRing<Ext> extensionPolynomialRing = (AbstractPolynomialRing<Ext>) AbstractPolynomialRing
				.getPolynomialRing(extension.extension(), numberOfVariables(), getComparator());
		Polynomial<Ext> extensionPolynomial = extensionPolynomialRing.getEmbedding(t, extension.embeddingMap());
		SquareFreeFactorizationResult<Ext> factors = extensionPolynomialRing
				.factorizeSquareFreePolynomial(extensionPolynomial);
		Set<Polynomial<Ext>> used = new TreeSet<>();
		List<Polynomial<Ext>> reduced = new ArrayList<>();
		MathMap<Ext, Ext> frobenious = new MathMap<>() {
			@Override
			public Ext evaluate(Ext t) {
				return extension.extension().power(t, getRing().getNumberOfElements());
			}
		};
		Ext unit = factors.unit;
		for (Polynomial<Ext> factor : factors.factors) {
			if (used.contains(factor)) {
				continue;
			}
			Polynomial<Ext> second = extensionPolynomialRing.getEmbedding(factor, frobenious);
			QuotientAndRemainderResult<Polynomial<Ext>> qr = extensionPolynomialRing.quotientAndRemainder(factor,
					second);
			if (qr.getRemainder().equals(extensionPolynomialRing.zero())) {
				Ext hilbert90 = extension.extension().hilbert90(qr.getQuotient().leadingCoefficient());
				unit = extension.extension().multiply(hilbert90, unit);
				reduced.add(extensionPolynomialRing.divideScalar(factor, hilbert90));
			} else {
				Polynomial<Ext> multiplied = extensionPolynomialRing.multiply(factor, second);
				Polynomial<Ext> multipliedSigma = extensionPolynomialRing.getEmbedding(multiplied, frobenious);
				Ext multiplier = extensionPolynomialRing.divideChecked(multiplied, multipliedSigma)
						.leadingCoefficient();
				Ext hilbert90 = extension.extension().hilbert90(multiplier);
				used.add(second);
				unit = extension.extension().multiply(hilbert90, unit);
				reduced.add(extensionPolynomialRing.divideScalar(multiplied, hilbert90));
			}
		}
		List<Polynomial<T>> result = new ArrayList<>();
		for (Polynomial<Ext> realFactor : reduced) {
			result.add(getEmbedding(realFactor, extension.retractionMap()));
		}
		SquareFreeFactorizationResult<T> r = new SquareFreeFactorizationResult<>();
		r.unit = extension.retractionMap().evaluate(unit);
		r.factors = result;
		return r;
	}

	@Override
	public FactorizationResult<Polynomial<Polynomial<T>>, Polynomial<T>> factorization(
			UnivariatePolynomial<Polynomial<T>> t) {
		PolynomialRing<T> addVar = AbstractPolynomialRing.getPolynomialRing(getRing(), numberOfVariables() + 1,
				getComparator());
		FactorizationResult<Polynomial<T>, Polynomial<T>> factors = addVar
				.uniqueFactorization(addVar.fromUnivariatePolynomial(t, numberOfVariables() + 1));
		SortedMap<Polynomial<Polynomial<T>>, Integer> result = new TreeMap<>();
		for (Polynomial<T> factor : factors.primeFactors()) {
			result.put(addVar.asUnivariatePolynomial(factor, numberOfVariables() + 1), factors.multiplicity(factor));
		}
		return new FactorizationResult<>(getEmbedding(factors.getUnit()), result);
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> squareFreeFactorization(Polynomial<T> t) {
		SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
		UnivariatePolynomialRing<Polynomial<T>> asYPolynomialRing = eliminateVariable().getUnivariatePolynomialRing();
		UnivariatePolynomial<Polynomial<T>> asYPolynomial = asUnivariatePolynomial(t, numberOfVariables());
		Polynomial<T> content = asYPolynomialRing.content(asYPolynomial);
		FactorizationResult<Polynomial<T>, T> contentFactors = eliminateVariable().squareFreeFactorization(content);
		for (Polynomial<T> factor : contentFactors.primeFactors()) {
			result.put(getEmbedding(factor), contentFactors.multiplicity(factor));
		}
		UnivariatePolynomial<Polynomial<T>> contentFree = asYPolynomialRing.contentFree(asYPolynomial);
		FactorizationResult<Polynomial<Polynomial<T>>, Polynomial<T>> squareFree = eliminateVariable()
				.getUnivariatePolynomialRing().squareFreeFactorization(contentFree);
		for (Polynomial<Polynomial<T>> factor : squareFree.primeFactors()) {
			result.put(fromUnivariatePolynomial(eliminateVariable().getUnivariatePolynomialRing().toUnivariate(factor),
					numberOfVariables()), squareFree.multiplicity(factor));
		}
		T unit = contentFactors.getUnit();
		unit = getRing().multiply(squareFree.getUnit().leadingCoefficient(), unit);
		return new FactorizationResult<>(unit, result);
	}

	@Override
	public FactorizationResult<Polynomial<T>, Polynomial<T>> uniqueFactorization(Polynomial<T> t) {
		if (t.degree() < 0) {
			throw new ArithmeticException("Can not factor zero!");
		}
		if (t.degree() == 0) {
			FactorizationResult<T, T> factors = getRing().uniqueFactorization(t.leadingCoefficient());
			SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
			for (T factor : factors.primeFactors()) {
				result.put(getEmbedding(factor), factors.multiplicity(factor));
			}
			return new FactorizationResult<>(getEmbedding(factors.getUnit()), result);
		}
		if (numberOfVariables() == 0) {
			throw new ArithmeticException("zero variable polynomial with degree > 0!");
		}
		if (numberOfVariables() == 1) {
			FactorizationResult<Polynomial<T>, T> factors = getRing()
					.factorization(getRing().getUnivariatePolynomialRing().toUnivariate(t));
			SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
			for (Polynomial<T> factor : factors.primeFactors()) {
				result.put(factor, factors.multiplicity(factor));
			}
			return new FactorizationResult<Polynomial<T>, Polynomial<T>>(getEmbedding(factors.getUnit()), result);
		}
		SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
		// UnivariatePolynomialRing<Polynomial<T>> asYPolynomialRing =
		// eliminateVariable().getUnivariatePolynomialRing();
		// UnivariatePolynomial<Polynomial<T>> asYPolynomial = asUnivariatePolynomial(t,
		// numberOfVariables());
		// Polynomial<T> content = asYPolynomialRing.content(asYPolynomial);
//		FactorizationResult<Polynomial<T>> contentFactors = eliminateVariable().uniqueFactorization(content);
//		for (Polynomial<T> factor : contentFactors.primeFactors()) {
//			result.put(getEmbedding(factor), contentFactors.multiplicity(factor));
//		}
		// UnivariatePolynomial<Polynomial<T>> contentFree =
		// asYPolynomialRing.contentFree(asYPolynomial);
		FactorizationResult<Polynomial<T>, T> squareFree = squareFreeFactorization(t);
		T unit = squareFree.getUnit();
		for (Polynomial<T> squareFreeFactor : squareFree.primeFactors()) {
			int multiplicity = squareFree.multiplicity(squareFreeFactor);
			int var = -1;
			for (int i = 0; i < numberOfVariables(); i++) {
				Polynomial<T> derivative = derivative(squareFreeFactor, i + 1);
				if (!derivative.equals(zero())) {
					var = i + 1;
					break;
				}
			}
			if (var == -1) {
				throw new ArithmeticException("Polynomial not square free after square free factorization!");
			}
			int[] map = new int[numberOfVariables()];
			map[0] = var - 1;
			for (int i = 1; i < numberOfVariables(); i++) {
				if (i != var - 1) {
					map[i] = i;
				} else {
					map[i] = 0;
				}
			}
			squareFreeFactor = getEmbedding(squareFreeFactor, map);
			SquareFreeFactorizationResult<T> factors = factorizeSquareFreePolynomial(squareFreeFactor);
			unit = getRing().multiply(getRing().power(factors.unit, multiplicity), unit);
			for (Polynomial<T> factor : factors.factors) {
				result.put(getEmbedding(factor, map), multiplicity);
			}
		}
		return new FactorizationResult<>(getEmbedding(unit), result);
//		if (t.degree() == 0) {
//			FactorizationResult<T> result = getRing().uniqueFactorization(t.leadingCoefficient());
//			SortedMap<Polynomial<T>, Integer> embedded = new TreeMap<>();
//			for (T p : result.primeFactors()) {
//				embedded.put(getEmbedding(p), result.multiplicity(p));
//			}
//			return new FactorizationResult<>(getEmbedding(result.getUnit()), embedded);
//		}
//		if (numberOfVariables() == 2 && getRing() instanceof Field<?>) {
//			Field<T> r = (Field<T>) getRing();
//			UnivariatePolynomialRing<T> uni = r.getUnivariatePolynomialRing();
//			List<Polynomial<T>> factors = new ArrayList<>();
//			if (!t.leadingCoefficient().equals(r.one())) {
//				factors.add(getEmbedding(t.leadingCoefficient()));
//				t = normalize(t);
//			}
//			int accuracy = 2 * t.degree() + 1;
//			FormalPowerSeries<T> ps = new FormalPowerSeries<T>(r, accuracy);
//			UnivariatePolynomialRing<PowerSeries<T>> psR = ps.getUnivariatePolynomialRing();
//			UnivariatePolynomial<Polynomial<T>> asUnivariate = asUnivariatePolynomial(t, 1);
//			UnivariatePolynomial<PowerSeries<T>> f = psR.getEmbedding(asUnivariate, new MathMap<>() {
//				@Override
//				public PowerSeries<T> evaluate(Polynomial<T> t) {
//					return ps.getEmbedding(uni.toUnivariate(t));
//				}
//			});
//
//			// PowerSeries<T> lead = f.leadingCoefficient();
//			f = psR.normalize(f);
//			FactorizationResult<Polynomial<PowerSeries<T>>> factorization = ps.factorization(f);
//			List<Polynomial<PowerSeries<T>>> liftedFactors = new ArrayList<>();
//			for (Polynomial<PowerSeries<T>> factor : factorization.primeFactors()) {
//				if (factor.degree() > 0) {
//					liftedFactors.add(factor);
//				}
//			}
//			Map<Integer, Polynomial<PowerSeries<T>>> powerSeriesFactors = new TreeMap<>();
//			for (int i = 0; i < liftedFactors.size(); i++) {
//				powerSeriesFactors.put(i, liftedFactors.get(i));
//			}
//			List<CombinedFactors<T>> powerSeriesCombinedFactors = new ArrayList<>();
//			powerSeriesCombinedFactors.add(new CombinedFactors<T>(Collections.emptySortedSet(), psR.one()));
//			List<Polynomial<T>> polynomialFactors = new ArrayList<>();
//			while (powerSeriesFactors.size() > 0) {
//				CheckCombinationsResult<T> result = checkCombinations(t, powerSeriesFactors, powerSeriesCombinedFactors,
//						ps);
//				t = result.cofactor;
//				polynomialFactors.addAll(result.factors);
//				for (int k : result.usedFactors) {
//					powerSeriesFactors.remove(k);
//				}
//				powerSeriesCombinedFactors = result.combined;
//			}
//			SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
//			for (Polynomial<T> factor : polynomialFactors) {
//				if (!result.containsKey(factor)) {
//					result.put(factor, 0);
//				}
//				result.put(factor, result.get(factor) + 1);
//			}
//			return new FactorizationResult<>(one(), result);
//		}
//		throw new ArithmeticException("Could not factorize polynomial!");
	}

	private static class CheckCombinationsResult<T extends Element<T>> {
		private List<Polynomial<T>> factors = new ArrayList<>();
		private Polynomial<T> cofactor;
		private Set<Integer> usedFactors = new TreeSet<>();
		private List<CombinedFactors<T>> combined = new ArrayList<>();
	}

	private static class CombinedFactors<T extends Element<T>> {
		private SortedSet<Integer> usedIndeces;
		private Polynomial<T> combined;

		public CombinedFactors(SortedSet<Integer> usedIndeces, Polynomial<T> combined) {
			this.usedIndeces = usedIndeces;
			this.combined = combined;
		}
	}

	private CheckCombinationsResult<T> checkCombinations(Polynomial<T> t,
			Map<Integer, Polynomial<T>> powerSeriesFactors, List<CombinedFactors<T>> powerSeriesCombinedFactors,
			int accuracy) {
		CheckCombinationsResult<T> result = new CheckCombinationsResult<>();
		result.cofactor = t;
		UnivariatePolynomial<T> leading = getRing().getUnivariatePolynomialRing()
				.toUnivariate(asUnivariatePolynomial(t, 1).leadingCoefficient());
		for (int i : powerSeriesFactors.keySet()) {
			Polynomial<T> powerSeriesFactor = powerSeriesFactors.get(i);
			for (CombinedFactors<T> powerSeriesCombinedFactor : powerSeriesCombinedFactors) {
				if (powerSeriesCombinedFactor.usedIndeces.size() == powerSeriesFactors.size()) {
					System.err.println(t);
					throw new ArithmeticException("No combinations left!");
				}
				if (powerSeriesCombinedFactor.usedIndeces.size() != 0
						&& powerSeriesCombinedFactor.usedIndeces.first() <= i) {
					continue;
				}
				SortedSet<Integer> indeces = new TreeSet<>();
				indeces.addAll(powerSeriesCombinedFactor.usedIndeces);
				indeces.add(i);
				Polynomial<T> newCombined = moduloVarPower(
						multiply(powerSeriesFactor, powerSeriesCombinedFactor.combined), 2, accuracy);
				CheckFactorResult<T> cfr = checkFactor(result.cofactor, leading, newCombined, accuracy);
				if (cfr.foundFactor) {
					result.factors.add(cfr.factor);
					result.cofactor = cfr.cofactor;
					result.usedFactors.addAll(indeces);
					leading = getRing().getUnivariatePolynomialRing()
							.toUnivariate(asUnivariatePolynomial(result.cofactor, 1).leadingCoefficient());
					break;
				} else {
					result.combined.add(new CombinedFactors<T>(indeces, newCombined));
				}
			}
		}
		return result;
	}

	private static class CheckFactorResult<T extends Element<T>> {
		private boolean foundFactor = false;
		private Polynomial<T> factor = null;
		private Polynomial<T> cofactor = null;
	}

	private CheckFactorResult<T> checkFactor(Polynomial<T> t, UnivariatePolynomial<T> leading,
			Polynomial<T> potentialFactor, int accuracy) {
		CheckFactorResult<T> result = new CheckFactorResult<>();
		UnivariatePolynomialRing<Polynomial<T>> ring = eliminateVariable().getUnivariatePolynomialRing();
		Polynomial<T> potentialLeading = getEmbedding(getRing().getUnivariatePolynomialRing().divide(leading,
				asUnivariatePolynomial(potentialFactor, 1).leadingCoefficient()), new int[] { 1 });
		potentialFactor = multiply(potentialLeading, potentialFactor);
		Polynomial<T> factor = moduloVarPower(potentialFactor, 2, accuracy);
		factor = fromUnivariatePolynomial(ring.contentFree(asUnivariatePolynomial(factor, 1)), 1);
		QuotientAndRemainderResult<Polynomial<T>> qr = quotientAndRemainder(t, factor);
		if (qr.getRemainder().equals(zero())) {
			result.foundFactor = true;
			result.factor = factor;
			result.cofactor = qr.getQuotient();
		}
		return result;
	}

	@Override
	public List<AffinePoint<T>> solve(List<Polynomial<T>> polynomials) {
		if (!(getRing() instanceof Field<?>)) {
			throw new UnsupportedOperationException();
		}
		Field<T> field = (Field<T>) getRing();
		if (polynomials.isEmpty()) {
			throw new InfinityException();
		}
		if (numberOfVariables() == 1) {
			Set<T> roots = new TreeSet<>();
			roots.addAll(field.roots(polynomials.get(0)).keySet());
			for (int i = 1; i < polynomials.size(); i++) {
				roots.retainAll(field.roots(polynomials.get(i)).keySet());
			}
			List<AffinePoint<T>> result = new ArrayList<>();
			for (T root : roots) {
				result.add(new AffinePoint<T>(field, Collections.singletonList(root)));
			}
			return result;
		}
		if (getComparator() != Monomial.REVLEX) {
			PolynomialRing<T> revlexRing = getPolynomialRing(getRing(), numberOfVariables(), Monomial.REVLEX);
			List<Polynomial<T>> embeded = new ArrayList<>();
			int[] map = new int[numberOfVariables()];
			for (int i = 0; i < numberOfVariables(); i++) {
				map[i] = i;
			}
			for (Polynomial<T> polynomial : polynomials) {
				embeded.add(revlexRing.getEmbedding(polynomial, map));
			}
			return revlexRing.solve(embeded);
		}
		Ideal<Polynomial<T>> eliminationIdeal = getIdeal(polynomials);
		List<Polynomial<T>> generators = eliminationIdeal.generators();
		int variable = 0;
		int highest = 0;
		List<AffinePoint<T>> prevPartialSolutions = new ArrayList<>();
		List<AffinePoint<T>> partialSolutions = new ArrayList<>();
		partialSolutions.add(new AffinePoint<>(field, Collections.emptyList()));
		int[] map = new int[numberOfVariables()];
		for (int i = 0; i < map.length; i++) {
			map[i] = -1;
		}
		for (int i = generators.size() - 1; i >= 0; i--) {
			Polynomial<T> p = generators.get(i);
			highest = highestVariable(p);
			if (variable > 0) {
				map[variable - 1] = -1;
			}
			map[highest - 1] = 0;
			variable = highest;
			prevPartialSolutions.clear();
			prevPartialSolutions.addAll(partialSolutions);
			partialSolutions.clear();
			for (AffinePoint<T> prevSolution : prevPartialSolutions) {
				List<T> values = new ArrayList<>();
				values.addAll(prevSolution.getCoords());
				for (int j = values.size(); j < numberOfVariables(); j++) {
					values.add(null);
				}
				Polynomial<T> evaluated = partiallyEvaluate(p, values);
				if (highest <= prevSolution.getCoords().size() && prevSolution.getCoords().get(highest - 1) != null) {
					if (evaluated.equals(zero())) {
						partialSolutions.add(prevSolution);
					}
					continue;
				}
				Set<T> roots;
				if (evaluated.equals(zero())) {
					roots = Collections.singleton(null);
				} else {
					Polynomial<T> univar = field.getUnivariatePolynomialRing()
							.normalize(field.getUnivariatePolynomialRing().getEmbedding(evaluated, map));
					roots = field.roots(univar).keySet();
				}
				for (T root : roots) {
					List<T> coords = new ArrayList<>();
					coords.addAll(prevSolution.getCoords());
					if (highest <= coords.size()) {
						if (coords.get(highest - 1) != null) {
							throw new ArithmeticException("Something went wrong");
						}
						coords.set(highest - 1, root);
					} else {
						coords.add(root);
					}
					partialSolutions.add(new AffinePoint<>(field, coords));
				}
			}
		}
		if (highest != numberOfVariables()) {
			throw new ArithmeticException("Too many degrees of freedom to solve!");
		}
		return partialSolutions;
	}

	private int highestVariable(Polynomial<T> t) {
		int highest = 0;
		for (Monomial m : t.monomials()) {
			if (t.coefficient(m).equals(getRing().zero())) {
				continue;
			}
			for (int i = numberOfVariables(); i >= 1; i--) {
				if (m.exponents()[i - 1] != 0) {
					highest = Math.max(highest, i);
				}
			}
		}
		return highest;
	}
}
