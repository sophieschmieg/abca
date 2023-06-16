package fields.polynomials;

import java.io.IOException;
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
import fields.helper.FieldEmbedding;
import fields.helper.FieldOfFractions;
import fields.helper.FieldOfFractions.Fraction;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
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
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.DifferentialForms.DifferentialForm;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import fields.polynomials.Monomial.EliminationOrder;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.ConCatMap;
import util.MiscAlgorithms;
import util.PeekableReader;
import util.SingletonSortedMap;
import util.Trie;
import util.Trie.Node;
import varieties.affine.AffinePoint;

public abstract class AbstractPolynomialRing<T extends Element<T>> extends AbstractAlgebra<T, Polynomial<T>>
		implements PolynomialRing<T> {
	private String[] variableNames;
	private Trie variableTrie;
	private Map<String, Integer> variableMap;
	private DifferentialForms<T> differentialForms;
	private Map<Integer, PolynomialRing<T>> eliminatedVariables;

	protected AbstractPolynomialRing(boolean generateUnivariatePolynomialRing) {
		super(generateUnivariatePolynomialRing);
	}

	protected AbstractPolynomialRing(boolean generateUnivariatePolynomialRing, String[] variableNames) {
		super(generateUnivariatePolynomialRing);
		setVariableNames(variableNames);
	}

	public static <T extends Element<T>> PolynomialRing<T> getPolynomialRing(Ring<T> ring,
			Comparator<Monomial> comparator, String[] variableNames) {
		PolynomialRing<T> result = getPolynomialRing(ring, variableNames.length, comparator);
		result.setVariableNames(variableNames);
		return result;
	}

	public static <T extends Element<T>> PolynomialRing<T> getPolynomialRing(Ring<T> ring, int numvars,
			Comparator<Monomial> comparator) {
		if (numvars == 1) {
			return ring.getUnivariatePolynomialRing();
		}
		return new MultivariatePolynomialRing<T>(ring, numvars, comparator);
	}

	@Override
	public String[] getVariableNames() {
		return variableNames;
	}

	@Override
	public void setVariableNames(String[] variableNames) {
		if (variableNames.length != numberOfVariables()) {
			throw new ArithmeticException("mismatched variable names");
		}
		this.variableNames = variableNames;
		variableMap = null;
		variableTrie = null;
	}

	@Override
	public final PolynomialRing<T> eliminateVariable(int variable) {
		if (eliminatedVariables == null) {
			eliminatedVariables = new TreeMap<>();
		}
		if (!eliminatedVariables.containsKey(variable)) {
			String[] leftOverVariableNames = new String[numberOfVariables() - 1];
			int i = 0;
			for (int j = 0; j < variableNames.length; j++) {
				if (j == variable - 1) {
					continue;
				}
				leftOverVariableNames[i] = variableNames[j];
				i++;
			}
			eliminatedVariables.put(variable,
					AbstractPolynomialRing.getPolynomialRing(getRing(), getComparator(), leftOverVariableNames));
		}
		return eliminatedVariables.get(variable);
	}

	@Override
	public final PolynomialRing<T> eliminateVariable() {
		return eliminateVariable(1);
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
	public Ideal<T> annihilator() {
		return getRing().getZeroIdeal();
	}

	@Override
	public Polynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values) {
		Polynomial<T> result = this.zero();
		for (Monomial m : t.monomials()) {
			Polynomial<T> value = getEmbedding(t.coefficient(m));
			for (int i = 0; i < t.numberOfVariables(); i++) {
				value = multiply(value, power(values.get(i), m.exponents()[i]));
			}
			result = add(result, value);
		}
		return result;
	}

	@Override
	public List<Vector<T>> getSyzygies() {
		return Collections.emptyList();
	}

	private <S extends Element<S>> FieldOfFractionsResult<Polynomial<T>, TExt<S>> fieldOfFractions(
			FieldOfFractionsResult<T, S> baseFractions) {
		PolynomialRing<S> polynomialRing = AbstractPolynomialRing.getPolynomialRing(baseFractions.getField(),
				numberOfVariables(), getComparator());
		TranscendentalFieldExtension<S> fieldOfFractions = new TranscendentalFieldExtension<>(baseFractions.getField(),
				polynomialRing);
		MathMap<TExt<S>, S> denominator = new MathMap<>() {

			@Override
			public S evaluate(TExt<S> t) {
				Polynomial<S> numerator = t.getNumerator();
				T denominator = getRing().one();
				for (Monomial m : numerator.monomials()) {
					denominator = getRing().lcm(denominator,
							baseFractions.getDenominator().evaluate(numerator.coefficient(m)));
				}
				Polynomial<S> denom = t.getDenominator();
				for (Monomial m : denom.monomials()) {
					denominator = getRing().lcm(denominator,
							baseFractions.getDenominator().evaluate(denom.coefficient(m)));
				}
				return baseFractions.getEmbedding().evaluate(denominator);
			}
		};
		return new FieldOfFractionsResult<>(this, fieldOfFractions, new MathMap<>() {

			@Override
			public TExt<S> evaluate(Polynomial<T> t) {
				return fieldOfFractions.getEmbedding(polynomialRing.getEmbedding(t, baseFractions.getEmbedding()));
			}
		}, new MathMap<>() {

			@Override
			public Polynomial<T> evaluate(TExt<S> t) {
				Polynomial<S> numerator = t.getNumerator();
				return getEmbedding(polynomialRing.multiply(denominator.evaluate(t), numerator),
						baseFractions.getAsInteger());
			}
		}, new MathMap<>() {

			@Override
			public Polynomial<T> evaluate(TExt<S> t) {
				Polynomial<S> denom = t.getDenominator();
				return getEmbedding(polynomialRing.multiply(denominator.evaluate(t), denom),
						baseFractions.getAsInteger());
			}
		}, new MathMap<>() {

			@Override
			public Polynomial<T> evaluate(TExt<S> t) {
				return getEmbedding(t.asInteger(), baseFractions.getAsInteger());
			}
		});
	}

	@Override
	public FieldOfFractionsResult<Polynomial<T>, ?> fieldOfFractions() {
		return fieldOfFractions(getRing().fieldOfFractions());
	}

	@Override
	public boolean isSquareFree(Polynomial<T> t) {
		return squareFreeFactorization(t).squareFree();
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
	public Polynomial<T> multiply(T t1, Polynomial<T> t2, Polynomial<T> t3) {
		return multiply(getEmbedding(t1), t2, t3);
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

	@Override
	public PolynomialIdeal<T> homogenizeIdeal(Ideal<Polynomial<T>> ideal) {
		List<Polynomial<T>> homogenousGenerators = new ArrayList<>();
		PolynomialRing<T> homogenousPolynomialRing = AbstractPolynomialRing.getPolynomialRing(getRing(),
				numberOfVariables() + 1, getComparator());
		for (Polynomial<T> generator : ideal.generators()) {
			homogenousGenerators.add(homogenousPolynomialRing
					.homogenize(homogenousPolynomialRing.getEmbedding(generator), numberOfVariables() + 1));
		}
		return homogenousPolynomialRing.getIdeal(homogenousGenerators);
	}
	
	@Override
	public Set<Polynomial<T>> orbit(Polynomial<T> t) {
		List<Integer> indeces = new ArrayList<>();
		for (int i = 0; i < numberOfVariables(); i++) {
			indeces.add(i);
		}
		Set<Polynomial<T>> result = new TreeSet<>();
		for (List<Integer> permutation : MiscAlgorithms.permutations(indeces) ) {
			int[] map = new int[numberOfVariables()];
			for (int i = 0; i < numberOfVariables(); i++) {
				map[i] = permutation.get(i);
			}
			result.add(getEmbedding(t, map));
		}
		return result;
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
	public List<PolynomialIdeal<T>> maximalPrimeIdealChain() {
		return maximalPrimeIdealChain(getNilRadical());
	}

	@Override
	public List<PolynomialIdeal<T>> maximalPrimeIdealChain(Ideal<Polynomial<T>> start) {
		PolynomialIdeal<T> ideal = start.isPrime() ? (PolynomialIdeal<T>) start
				: primaryDecomposition(start).getRadicals().get(0);
		Ideal<T> intersectToRing = ideal.intersectToRing();
		List<? extends Ideal<T>> ringChain = getRing().maximalPrimeIdealChain(intersectToRing);
		List<PolynomialIdeal<T>> result = new ArrayList<>();
		for (Ideal<T> ringIdeal : ringChain) {
			ideal = add(ideal, getEmbeddingOfBaseIdeal(ringIdeal));
			result.add(ideal);
		}
		Set<Integer> boundVariables = ideal.divideOut().boundVariables();
		for (int i = 0; i < numberOfVariables(); i++) {
			if (!boundVariables.contains(i + 1)) {
				ideal = add(ideal, getIdeal(Collections.singletonList(getVar(i + 1))));
				result.add(ideal);
			}
		}
		return result;
	}

	@Override
	public List<PolynomialIdeal<T>> maximalPrimeIdealChain(Ideal<Polynomial<T>> start, Ideal<Polynomial<T>> end) {
		if (!end.contains(start) || !end.isPrime()) {
			throw new ArithmeticException("Preconditions not met!");
		}
		PolynomialIdeal<T> ideal = primaryDecomposition(start).getRadicals().get(0);
		Ideal<T> intersectToRing = ideal.intersectToRing();
		List<? extends Ideal<T>> ringChain = getRing().maximalPrimeIdealChain(intersectToRing);
		List<PolynomialIdeal<T>> result = new ArrayList<>();
		for (Ideal<T> ringIdeal : ringChain) {
			ideal = add(ideal, getEmbeddingOfBaseIdeal(ringIdeal));
			result.add(ideal);
		}
		for (Polynomial<T> generator : end.generators()) {
			if (!ideal.contains(generator)) {
				ideal = add(ideal, getIdeal(Collections.singletonList(generator)));
				result.add(ideal);
			}
		}
		return result;
	}

	@Override
	public QuotientAndRemainderResult<Polynomial<T>> quotientAndRemainder(Polynomial<T> dividend,
			Polynomial<T> divisor) {
		GeneralQuotientAndRemainderResult<T> gqr = generalQuotientAndRemainder(dividend,
				Collections.singletonList(divisor));
		new QuotientAndRemainderResult<>(gqr.getQuotients().get(0), gqr.getRemainder());
		QuotientAndRemainderResult<Polynomial<T>> result = new QuotientAndRemainderResult<>(gqr.getQuotients().get(0),
				gqr.getRemainder());
		return result;
	}

	public GeneralQuotientAndRemainderResult<T> generalQuotientAndRemainder(Polynomial<T> polynomial,
			List<Polynomial<T>> basis) {
		Ring<T> ring = getRing();
		Comparator<Monomial> comparator = getComparator();
		List<SortedMap<Monomial, T>> quotient = new ArrayList<>();
		SortedMap<Monomial, T> remainder = new TreeMap<>();
		SortedMap<Monomial, T> p = new TreeMap<>();
		for (Monomial m : polynomial.monomials()) {
			p.put(m, polynomial.coefficient(m));
		}
		for (int i = 0; i < basis.size(); i++) {
			quotient.add(new TreeMap<>(comparator));
		}
		while (true) {
			if (p.isEmpty()) {
				List<Polynomial<T>> list = new ArrayList<>();
				for (SortedMap<Monomial, T> q : quotient) {
					list.add(this.getPolynomial(q));
				}
				return new GeneralQuotientAndRemainderResult<>(this.getPolynomial(remainder), list);
			}
			Monomial leadingmonomial = p.lastKey();
			T leadingcoefficient = p.get(leadingmonomial);
			boolean success = false;
			boolean totalFailure = true;
			for (int i = 0; i < basis.size(); i++) {
				Polynomial<T> base = basis.get(i);
				if (comparator.compare(leadingmonomial, base.leadingMonomial()) >= 0
						&& !base.leadingCoefficient().equals(ring.zero())) {
					totalFailure = false;
				}
				Monomial div = leadingmonomial.divide(base.leadingMonomial());
				QuotientAndRemainderResult<T> qr = ring.quotientAndRemainder(leadingcoefficient,
						base.leadingCoefficient());
				if (div == null || qr.getQuotient().equals(ring.zero())) {
					continue;
				}
				T factor = qr.getQuotient();
				quotient.get(i).put(div, factor);
				for (Monomial m : base.monomials()) {
					T subtrahend = ring.multiply(factor, base.coefficient(m));
					Monomial multiplied = m.multiply(div);
					T newValue;
					if (p.containsKey(multiplied)) {
						newValue = ring.subtract(p.get(multiplied), subtrahend);
					} else {
						newValue = ring.negative(subtrahend);
					}
					if (newValue.equals(ring.zero())) {
						p.remove(multiplied);
					} else {
						p.put(multiplied, newValue);
					}
				}
				success = true;
				break;
			}
			if (totalFailure) {
				remainder.putAll(p);
				p.clear();
			} else if (!success) {
				remainder.put(leadingmonomial, leadingcoefficient);
				p.remove(leadingmonomial);
			}
		}
	}

	public static class ReduceAndExpressResult<T extends Element<T>> implements Comparable<ReduceAndExpressResult<T>> {
		private Polynomial<T> reduced;
		private List<Polynomial<T>> express;
		private boolean modified;

		private ReduceAndExpressResult(Polynomial<T> reduced, List<Polynomial<T>> express, boolean modified) {
			this.reduced = reduced;
			this.express = express;
			this.modified = modified;
		}

		public Polynomial<T> getReduced() {
			return reduced;
		}

		public List<Polynomial<T>> getExpression() {
			return express;
		}

		public boolean isModified() {
			return modified;
		}

		@Override
		public int compareTo(ReduceAndExpressResult<T> o) {
			return reduced.compareTo(o.reduced);
		}
	}

	private ReduceAndExpressResult<T> reduceAndExpress(Polynomial<T> polynomial,
			List<Polynomial<T>> polynomialExpression, List<ReduceAndExpressResult<T>> basis,
			boolean computeExpressionsAndSyzygies) {
		List<Polynomial<T>> extracted = new ArrayList<>();
		List<List<Polynomial<T>>> expressions = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			extracted.add(basis.get(i).getReduced());
			expressions.add(basis.get(i).getExpression());
		}
		return reduceAndExpress(polynomial, polynomialExpression, extracted, expressions,
				computeExpressionsAndSyzygies);
	}

	public ReduceAndExpressResult<T> reduceAndExpress(Polynomial<T> polynomial,
			List<Polynomial<T>> polynomialExpression, List<Polynomial<T>> basis,
			List<List<Polynomial<T>>> basisExpressions, boolean computeExpressionsAndSyzygies) {
		if (basis.isEmpty()) {
			Polynomial<T> unit = inverse(projectToUnit(polynomial));
			if (unit.equals(one())) {
				return new ReduceAndExpressResult<>(polynomial, polynomialExpression, false);
			}
			List<Polynomial<T>> expression = new ArrayList<>();
			if (computeExpressionsAndSyzygies) {
				for (Polynomial<T> coefficient : polynomialExpression) {
					expression.add(multiply(unit, coefficient));
				}
			}
			return new ReduceAndExpressResult<>(multiply(unit, polynomial), expression, true);
		}
		GeneralQuotientAndRemainderResult<T> gqr = generalQuotientAndRemainder(polynomial, basis);
		// f = a_1 b_1 + ... a_n_b_n + r' = a_1 (c_11 g_1 + c_12 g2 + ... c_1m g_m) +
		// ... + r'
		// f = d_1 g_1 + ... + d_m g_m
		// r' = ur
		// r = u^-1((d_1-a_1c11-a_2c21-...)g_1 + ...)
		Polynomial<T> unit = inverse(projectToUnit(gqr.getRemainder()));
		boolean modified = !(unit.equals(one()) && polynomial.equals(gqr.getRemainder()));
		List<Polynomial<T>> expression = new ArrayList<>();
		if (computeExpressionsAndSyzygies) {
			for (int j = 0; j < polynomialExpression.size(); j++) {
				expression.add(multiply(unit, polynomialExpression.get(j)));
			}
			for (int i = 0; i < basis.size(); i++) {
				Polynomial<T> quotient = gqr.getQuotients().get(i);
				List<Polynomial<T>> basisExpression = basisExpressions.get(i);
				if (computeExpressionsAndSyzygies) {
					for (int j = 0; j < expression.size(); j++) {
						expression.set(j,
								subtract(expression.get(j), multiply(unit, quotient, basisExpression.get(j))));
					}
				}
			}
		}
		return new ReduceAndExpressResult<>(multiply(unit, gqr.getRemainder()), expression, modified);
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
	public Vector<T> evaluate(Vector<Polynomial<T>> t, Vector<T> ts) {
		return Vector.mapVector(new MathMap<Polynomial<T>, T>() {
			@Override
			public T evaluate(Polynomial<T> t) {
				return AbstractPolynomialRing.this.evaluate(t, ts);
			}
		}, t);
	}

	@Override
	public Matrix<T> evaluate(Matrix<Polynomial<T>> t, Vector<T> ts) {
		return Matrix.mapMatrix(new MathMap<Polynomial<T>, T>() {
			@Override
			public T evaluate(Polynomial<T> t) {
				return AbstractPolynomialRing.this.evaluate(t, ts);
			}
		}, t);
	}

	@Override
	public Vector<Polynomial<T>> gradient(Polynomial<T> t) {
		return differentialForms().asGradedVector(totalDerivative(t), 1);
	}

	@Override
	public Matrix<Polynomial<T>> jacobianMatrix(Vector<Polynomial<T>> t) {
		List<Vector<Polynomial<T>>> rows = new ArrayList<>();
		for (int i = 0; i < t.dimension(); i++) {
			rows.add(gradient(t.get(i + 1)));
		}
		return Matrix.fromRows(rows);
	}

	@Override
	public DifferentialForm<T> totalDerivative(Polynomial<T> t) {
		return differentialForms().derivative(t);
	}

	@Override
	public DifferentialForms<T> differentialForms() {
		if (differentialForms == null) {
			differentialForms = new DifferentialForms<>(this);
		}
		return differentialForms;
	}

	@Override
	public final MultivariatePolynomialRing<T> addVariableWithElimination(int shift) {
		if (shift < 0)
			throw new ArithmeticException("Cannot count");
		return new MultivariatePolynomialRing<T>(getRing(), numberOfVariables() + shift,
				new EliminationOrder(Monomial.GREVLEX, getComparator(), shift));
	}

	@Override
	public PolynomialIdeal<T> getIdeal(List<Polynomial<T>> generators) {
		if (getRing() instanceof PolynomialRing<?>) {
			return getIdealWithTransforms(generators).getIdeal();
		}
		GroebnerBasis<T> basis = buchberger(generators, false);
		return new PolynomialIdeal<T>(this, basis.getBasis());
	}

	@Override
	public Polynomial<T> flattenPolynomial(Polynomial<Polynomial<T>> t) {
		Polynomial<T> result = zero();
		for (Monomial m : t.monomials()) {
			Polynomial<T> coefficient = getEmbedding(t.coefficient(m));
			int[] exponents = new int[numberOfVariables()];
			for (int i = 0; i < m.exponents().length; i++) {
				exponents[numberOfVariables() - m.exponents().length + i] = m.exponents()[i];
			}
			Monomial multiplier = getMonomial(exponents);
			result = add(result, multiply(getRing().one(), multiplier, coefficient));
		}
		return result;
	}

	@Override
	public Polynomial<Polynomial<T>> unflattenPolynomial(Polynomial<T> t, PolynomialRing<Polynomial<T>> polynomialRing,
			PolynomialRing<T> baseRing) {
		if (numberOfVariables() != baseRing.numberOfVariables() + polynomialRing.numberOfVariables()) {
			throw new ArithmeticException("Variables do not add up");
		}
		Polynomial<Polynomial<T>> result = polynomialRing.zero();
		for (Monomial m : t.monomials()) {
			int[] baseExponents = Arrays.copyOf(m.exponents(), baseRing.numberOfVariables());
			Monomial base = baseRing.getMonomial(baseExponents);
			Polynomial<T> baseCoeff = baseRing.getPolynomial(Collections.singletonMap(base, t.coefficient(m)));
			int[] exponents = new int[polynomialRing.numberOfVariables()];
			for (int i = baseRing.numberOfVariables(); i < numberOfVariables(); i++) {
				exponents[i - baseRing.numberOfVariables()] = m.exponents()[i];
			}
			result = polynomialRing.add(result,
					polynomialRing.multiply(baseCoeff, polynomialRing.getMonomial(exponents), polynomialRing.one()));
		}
		return result;
	}

	@SuppressWarnings("unchecked")
	private <S extends Element<S>> IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransformsOverPolynomialRing(
			List<Polynomial<T>> generators, PolynomialRing<S> ring) {
		PolynomialRing<S> combinedRing = AbstractPolynomialRing.getPolynomialRing(ring.getRing(),
				ring.numberOfVariables() + numberOfVariables(),
				new Monomial.EliminationOrder(ring.getComparator(), getComparator(), ring.numberOfVariables()));
		List<Polynomial<S>> combinedGenerators = new ArrayList<>();
		for (Polynomial<T> generator : generators) {
			combinedGenerators.add(combinedRing.flattenPolynomial((Polynomial<Polynomial<S>>) generator));
		}
		IdealResult<Polynomial<S>, PolynomialIdeal<S>> ideal = combinedRing.getIdealWithTransforms(combinedGenerators);
		List<Polynomial<T>> basis = new ArrayList<>();
		for (Polynomial<S> generator : ideal.getIdeal().generators()) {
			basis.add((Polynomial<T>) combinedRing.unflattenPolynomial(generator, (PolynomialRing<Polynomial<S>>) this,
					ring));
		}
		List<List<Polynomial<T>>> transforms = new ArrayList<>();
		for (List<Polynomial<S>> transform : ideal.getGeneratorExpressions()) {
			List<Polynomial<T>> transformList = new ArrayList<>();
			for (Polynomial<S> t : transform) {
				transformList.add((Polynomial<T>) combinedRing.unflattenPolynomial(t,
						(PolynomialRing<Polynomial<S>>) this, ring));
			}
			transforms.add(transformList);
		}
		List<Vector<Polynomial<T>>> syzygies = new ArrayList<>();
		for (Vector<Polynomial<S>> syzygy : ideal.getSyzygies()) {
			List<Polynomial<T>> syzygyList = new ArrayList<>();
			for (Polynomial<S> t : syzygy.asList()) {
				syzygyList.add((Polynomial<T>) combinedRing.unflattenPolynomial(t, (PolynomialRing<Polynomial<S>>) this,
						ring));
			}
			syzygies.add(new Vector<>(syzygyList));
		}
		return new IdealResult<>(transforms, generators, new PolynomialIdeal<>(this, basis), syzygies);
	}

	@Override
	public IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransforms(List<Polynomial<T>> generators) {
		if (getRing() instanceof PolynomialRing<?>) {
			return getIdealWithTransformsOverPolynomialRing(generators, (PolynomialRing<?>) getRing());
		}
		GroebnerBasis<T> basis = buchberger(generators, true);
		return new IdealResult<>(basis.getExpression(), generators, new PolynomialIdeal<T>(this, basis.getBasis()),
				basis.getSyzygies());
	}

	@Override
	public PolynomialIdeal<T> getIdeal(@SuppressWarnings("unchecked") Polynomial<T>... generators) {
		return getIdeal(Arrays.asList(generators));
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
	public Ideal<Polynomial<T>> getNilRadical() {
		return getEmbeddingOfBaseIdeal(getRing().getNilRadical());
	}

	@Override
	public PolynomialIdeal<T> getUnitIdeal() {
		return getIdeal(Collections.singletonList(one()));
	}

	@Override
	public FactorizationResult<Ideal<Polynomial<T>>, Ideal<Polynomial<T>>> idealFactorization(Ideal<Polynomial<T>> t) {
		if (numberOfVariables() == 0) {
			List<T> generators = new ArrayList<>();
			for (Polynomial<T> generator : t.generators()) {
				generators.add(generator.leadingCoefficient());
			}
			Ideal<T> ideal = getRing().getIdeal(generators);
			FactorizationResult<Ideal<T>, Ideal<T>> idealFactors = getRing().idealFactorization(ideal);
			SortedMap<Ideal<Polynomial<T>>, Integer> result = new TreeMap<>();
			for (Ideal<T> primeFactor : idealFactors.primeFactors()) {
				List<Polynomial<T>> factorGenerators = new ArrayList<>();
				for (T generator : primeFactor.generators()) {
					factorGenerators.add(getEmbedding(generator));
				}
				Ideal<Polynomial<T>> asPolynomialIdeal = getIdeal(factorGenerators);
				result.put(asPolynomialIdeal, idealFactors.multiplicity(primeFactor));
			}
			return new FactorizationResult<>(getUnitIdeal(), result);
		}
		if (!getRing().isIntegral() || getRing().krullDimension() != 0 || numberOfVariables() != 1) {
			return super.idealFactorization(t);
		}
		throw new ArithmeticException("Not a Dedekind ring!");
	}

	@Override
	public PolynomialIdeal<T> getEmbeddingOfBaseIdeal(Ideal<T> t) {
		List<Polynomial<T>> generators = new ArrayList<>();
		for (T generator : t.generators()) {
			generators.add(getEmbedding(generator));
		}
		return getIdeal(generators);
	}

	@Override
	public PolynomialIdeal<T> getEmbedding(Ideal<Polynomial<T>> t) {
		PolynomialIdeal<T> ideal = (PolynomialIdeal<T>) t;
		if (ideal.dimension() != -3/* 0 */) {
			List<Polynomial<T>> generators = new ArrayList<>();
			for (Polynomial<T> generator : t.generators()) {
				generators.add(getEmbedding(generator));
			}
			return getIdeal(generators);
		}
		PolynomialRing<T> ring = ideal.getRing();
		return null;
	}

	private List<Monomial> staircase(PolynomialIdeal<T> ideal) {
		List<Monomial> result = new ArrayList<>();
		CoordinateRing<T> cr = ideal.divideOut();
		Integers z = Integers.z();
		List<Monomial> leadingMonomials = new ArrayList<>();
		for (Polynomial<T> generator : ideal.generators()) {
			leadingMonomials.add(generator.leadingMonomial());
		}
		for (int i = 0; !cr.hilbertFunction(i).equals(z.zero()); i++) {
			for (Monomial m : monomials(i)) {
				boolean divisible = false;
				for (Monomial leading : leadingMonomials) {
					if (m.divide(leading) != null) {
						divisible = true;
						break;
					}
				}
				if (!divisible) {
					result.add(m);
				}
			}
		}
		return result;
	}

	private List<Monomial> monomials(int degree) {
		return monomials(numberOfVariables(), degree);
	}

	private List<Monomial> monomials(int variable, int degree) {
		int[] exp = new int[numberOfVariables()];
		if (variable == 0 && degree == 0) {
			return Collections.singletonList(getMonomial(exp));
		}
		if (variable == 0) {
			return Collections.emptyList();
		}
		if (variable == 1) {
			exp[0] = degree;
			return Collections.singletonList(getMonomial(exp));
		}
		List<Monomial> result = new ArrayList<>();
		for (int i = 0; i <= degree; i++) {
			exp[variable - 1] = i;
			Monomial m1 = getMonomial(exp);
			for (Monomial m2 : monomials(variable - 1, degree - i)) {
				result.add(m1.multiply(m2));
			}
		}
		return result;
	}

	@Override
	public PolynomialIdeal<T> getEmbedding(Ideal<Polynomial<T>> t, int[] map) {
		List<Polynomial<T>> generators = new ArrayList<>();
		for (Polynomial<T> generator : t.generators()) {
			generators.add(getEmbedding(generator, map));
		}
		return getIdeal(generators);
	}

	@Override
	public <S extends Element<S>> PolynomialIdeal<T> getEmbedding(Ideal<Polynomial<S>> t, MathMap<S, T> map) {
		List<Polynomial<T>> generators = new ArrayList<>();
		for (Polynomial<S> generator : t.generators()) {
			generators.add(getEmbedding(generator, map));
		}
		return getIdeal(generators);
	}

	private static <T extends Element<T>> Polynomial<Polynomial<T>> asPolynomialOverUnivariate(Polynomial<T> polynomial,
			UnivariatePolynomialRing<T> univariate, PolynomialRing<Polynomial<T>> polynomialRing) {
		Polynomial<Polynomial<T>> result = polynomialRing.zero();
		for (Monomial m : polynomial.monomials()) {
			int[] exp = Arrays.copyOf(m.exponents(), polynomialRing.numberOfVariables());
			Monomial newM = polynomialRing.getMonomial(exp);
			UnivariatePolynomial<T> univariateCoeff = univariate.multiply(polynomial.coefficient(m),
					univariate.getVarPower(m.exponents()[polynomialRing.numberOfVariables()]));
			Polynomial<Polynomial<T>> coeff = polynomialRing
					.getPolynomial(Collections.singletonMap(newM, univariateCoeff));
			result = polynomialRing.add(result, coeff);
		}
		return result;
	}

	private static <T extends Element<T>> Polynomial<T> fromPolynomialOverUnivariate(
			Polynomial<Polynomial<T>> polynomial, UnivariatePolynomialRing<T> univariate,
			PolynomialRing<T> polynomialRing) {
		Polynomial<T> result = polynomialRing.zero();
		for (Monomial m : polynomial.monomials()) {
			int[] exp = Arrays.copyOf(m.exponents(), polynomialRing.numberOfVariables());
			UnivariatePolynomial<T> coeff = univariate.toUnivariate(polynomial.coefficient(m));
			Map<Monomial, T> coeffs = new TreeMap<>();
			for (int i = 0; i <= coeff.degree(); i++) {
				exp[polynomialRing.numberOfVariables() - 1] = i;
				Monomial newM = polynomialRing.getMonomial(exp);
				coeffs.put(newM, coeff.univariateCoefficient(i));
			}
			result = polynomialRing.add(result, polynomialRing.getPolynomial(coeffs));
		}
		return result;
	}

	private <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, Ext extends FieldExtension<B, E, Ext>> PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionZeroDimensional(
			PolynomialIdeal<T> ideal, Ext ringByMax, MathMap<T, E> embedding, MathMap<E, T> inverse, Ideal<T> maximal) {
		if (ideal.dimension() != 0) {
			throw new ArithmeticException("Not a zero dimensional ideal!");
		}
		if (numberOfVariables() == 0) {
			List<Polynomial<T>> generators = new ArrayList<>();
			for (T generator : maximal.generators()) {
				generators.add(getEmbedding(generator));
			}
			return new PrimaryDecompositionResult<>(Collections.singletonList(ideal),
					Collections.singletonList(getIdeal(generators)));
		}
		if (getComparator() != Monomial.LEX) {
			AbstractPolynomialRing<T> revlexRing = (AbstractPolynomialRing<T>) AbstractPolynomialRing
					.getPolynomialRing(getRing(), numberOfVariables(), Monomial.LEX);
			PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> decomposition = revlexRing
					.primaryDecompositionZeroDimensional(revlexRing.getEmbedding(ideal), ringByMax, embedding, inverse,
							maximal);
			List<PolynomialIdeal<T>> primaries = new ArrayList<>();
			List<PolynomialIdeal<T>> radicals = new ArrayList<>();
			for (PolynomialIdeal<T> primary : decomposition.getPrimaryIdeals()) {
				primaries.add(getEmbedding(primary));
			}
			for (PolynomialIdeal<T> radical : decomposition.getRadicals()) {
				radicals.add(getEmbedding(radical));
			}
			return new PrimaryDecompositionResult<>(primaries, radicals);
		}
		List<Polynomial<T>> intersectToUnivariate = new ArrayList<>();
		generatorLoop: for (Polynomial<T> generator : ideal.generators()) {
			for (int i = 1; i < numberOfVariables(); i++) {
				if (generator.degree(i) != 0) {
					continue generatorLoop;
				}
			}
			intersectToUnivariate.add(generator);
		}
		PolynomialIdeal<T> intersectToUnivariateIdeal = getIdeal(intersectToUnivariate);
		Polynomial<T> maximalDegree = null;
		for (Polynomial<T> generator : intersectToUnivariateIdeal.generators()) {
			if (maximalDegree == null || generator.degree() > maximalDegree.degree()) {
				maximalDegree = generator;
			}
		}
		UnivariatePolynomial<E> reduced = ringByMax.getUnivariatePolynomialRing().getEmbedding(
				getRing().getUnivariatePolynomialRing().getEmbedding(maximalDegree, new int[numberOfVariables()]),
				embedding);
		UnivariatePolynomialRing<T> univariateRing = getRing().getUnivariatePolynomialRing();
		AbstractPolynomialRing<Polynomial<T>> newRing = (AbstractPolynomialRing<Polynomial<T>>) AbstractPolynomialRing
				.getPolynomialRing(univariateRing, numberOfVariables() - 1, Monomial.REVLEX);

		FactorizationResult<Polynomial<E>, E> factorization = ringByMax.factorization(reduced);
		Polynomial<T> product = one();
		for (Polynomial<E> factor : factorization.primeFactors()) {
			Polynomial<T> polynomial = getEmbedding(
					getRing().getUnivariatePolynomialRing().getEmbedding(factor, inverse),
					new int[] { numberOfVariables() - 1 });
			product = multiply(product, power(polynomial, factorization.multiplicity(factor)));
		}
		int exponent = 1;
		Polynomial<T> power = product;
		while (true) {
			if (intersectToUnivariateIdeal.contains(power)) {
				break;
			}
			exponent++;
			power = multiply(power, product);
		}
		List<PolynomialIdeal<T>> result = new ArrayList<>();
		List<PolynomialIdeal<T>> radicals = new ArrayList<>();
		for (Polynomial<E> factor : factorization.primeFactors()) {
			UnivariatePolynomial<T> univariate = getRing().getUnivariatePolynomialRing().getEmbedding(factor, inverse);
			Polynomial<T> polynomial = getEmbedding(univariate, new int[] { numberOfVariables() - 1 });
			int multiplicity = exponent * factorization.multiplicity(factor);
			List<Polynomial<T>> maximalGenerators = new ArrayList<>();
			maximalGenerators.add(univariate);
			for (T generator : maximal.generators()) {
				maximalGenerators.add(getRing().getUnivariatePolynomialRing().getEmbedding(generator));
			}
			PolynomialIdeal<T> newMaximal = getRing().getUnivariatePolynomialRing().getIdeal(maximalGenerators);
			FieldEmbedding<B, E, Ext> newExtension = ringByMax
					.getEmbeddedExtension(ringByMax.getUnivariatePolynomialRing().toUnivariate(factor));
			MathMap<Polynomial<T>, E> newEmbedding = new MathMap<>() {
				@Override
				public E evaluate(Polynomial<T> t) {
					return newExtension
							.fromPolynomial(ringByMax.getUnivariatePolynomialRing().getEmbedding(t, embedding));
				}
			};
			MathMap<E, Polynomial<T>> newInverse = new MathMap<>() {
				@Override
				public Polynomial<T> evaluate(E t) {
					return getRing().getUnivariatePolynomialRing().getEmbedding(newExtension.asPolynomial(t), inverse);
				}
			};
			List<Polynomial<Polynomial<T>>> newGenerators = new ArrayList<>();
			newGenerators.add(asPolynomialOverUnivariate(power(polynomial, multiplicity), univariateRing, newRing));
			for (Polynomial<T> generator : ideal.generators()) {
				newGenerators.add(asPolynomialOverUnivariate(generator, univariateRing, newRing));
			}
			PolynomialIdeal<Polynomial<T>> newIdeal = newRing.getIdeal(newGenerators);
			PrimaryDecompositionResult<Polynomial<Polynomial<T>>, PolynomialIdeal<Polynomial<T>>> decomposition = newRing
					.primaryDecompositionZeroDimensional(newIdeal, newExtension.getField(), newEmbedding, newInverse,
							newMaximal);
			for (PolynomialIdeal<Polynomial<T>> primary : decomposition.getPrimaryIdeals()) {
				List<Polynomial<T>> generators = new ArrayList<>();
				for (Polynomial<Polynomial<T>> generator : primary.generators()) {
					generators.add(fromPolynomialOverUnivariate(generator, univariateRing, this));
				}
				result.add(getIdeal(generators));
			}
			for (PolynomialIdeal<Polynomial<T>> radical : decomposition.getRadicals()) {
				List<Polynomial<T>> generators = new ArrayList<>();
				for (Polynomial<Polynomial<T>> generator : radical.generators()) {
					generators.add(fromPolynomialOverUnivariate(generator, univariateRing, this));
				}
				radicals.add(getIdeal(generators));
			}
		}
		return new PrimaryDecompositionResult<>(result, radicals);
	}

	private <S extends Element<S>, B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, Ext extends FieldExtension<B, E, Ext>, R extends Ring<T>, I extends Ideal<T>, F extends Field<S>> PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionZeroDimensional(
			PolynomialIdeal<T> t, Ideal<T> maximalIdeal, ModuloMaximalIdealResult<T, S, R, I, F> mod,
			Extension<S, B, E, Ext> extension) {
		return primaryDecompositionZeroDimensional(t, extension.extension(),
				new ConCatMap<>(mod.getReduction(), extension.embeddingMap()),
				new ConCatMap<>(extension.retractionMap(), mod.getLift()), maximalIdeal);
	}

	private <S extends Element<S>, R extends Ring<T>, I extends Ideal<T>, F extends Field<S>> PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionZeroDimensional(
			PolynomialIdeal<T> t, Ideal<T> maximalIdeal, ModuloMaximalIdealResult<T, S, R, I, F> mod) {
		return primaryDecompositionZeroDimensional(t, maximalIdeal, mod,
				mod.getField().getExtension(mod.getField().getUnivariatePolynomialRing().getVar()));
	}

	private PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionZeroDimensional(
			PolynomialIdeal<T> t, Ideal<T> maximalIdeal) {
		return primaryDecompositionZeroDimensional(t, maximalIdeal, getRing().moduloMaximalIdeal(maximalIdeal));
	}

	private T findDenominator(PolynomialIdeal<T> ideal, T primeIdealGenerator) {
		Ring<T> ring = getRing();
		T s = ring.one();
		if (primeIdealGenerator.equals(ring.zero())) {
			for (Polynomial<T> generator : ideal.generators()) {
				s = ring.multiply(s, generator.leadingCoefficient());
			}
			return s;
		}
		for (Polynomial<T> generator : ideal.generators()) {
			T lc = generator.leadingCoefficient();
			while (ring.isDivisible(lc, primeIdealGenerator)) {
				lc = ring.divideChecked(lc, primeIdealGenerator);
			}
			s = ring.multiply(s, lc);
		}
		return s;
	}

	private <S extends Element<S>> PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionZeroDimensionalFieldOfFractions(
			PolynomialIdeal<T> t, FieldOfFractionsResult<T, S> baseFractions) {
		PolynomialRing<S> polynomialRing = AbstractPolynomialRing.getPolynomialRing(baseFractions.getField(),
				numberOfVariables(), getComparator());
		PolynomialIdeal<S> ideal = polynomialRing.getEmbedding(t, baseFractions.getEmbedding());
		PrimaryDecompositionResult<Polynomial<S>, PolynomialIdeal<S>> decomposition = polynomialRing
				.primaryDecomposition(ideal);
		List<PolynomialIdeal<T>> primaries = new ArrayList<>();
		List<PolynomialIdeal<T>> radicals = new ArrayList<>();
		for (PolynomialIdeal<S> primary : decomposition.getPrimaryIdeals()) {
			List<Polynomial<T>> generators = new ArrayList<>();
			for (Polynomial<S> generator : primary.generators()) {
				T denominator = getRing().one();
				for (Monomial m : generator.monomials()) {
					denominator = getRing().lcm(denominator,
							baseFractions.getDenominator().evaluate(generator.coefficient(m)));
				}
				generators.add(contentFree(getEmbedding(
						polynomialRing.multiply(baseFractions.getEmbedding().evaluate(denominator), generator),
						baseFractions.getAsInteger())));
			}
			primaries.add(getIdeal(generators));
		}
		for (PolynomialIdeal<S> radical : decomposition.getRadicals()) {
			List<Polynomial<T>> generators = new ArrayList<>();
			for (Polynomial<S> generator : radical.generators()) {
				T denominator = getRing().one();
				for (Monomial m : generator.monomials()) {
					denominator = getRing().lcm(denominator,
							baseFractions.getDenominator().evaluate(generator.coefficient(m)));
				}
				generators.add(contentFree(getEmbedding(
						polynomialRing.multiply(baseFractions.getEmbedding().evaluate(denominator), generator),
						baseFractions.getAsInteger())));
			}
			radicals.add(getIdeal(generators));
		}
		return new PrimaryDecompositionResult<>(primaries, radicals);
	}

	private <S extends Element<S>, U extends Element<U>, R extends Element<R>> PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionZeroDimensionalLocalized(
			PolynomialIdeal<T> t, T primeIdealGenerator, LocalizeResult<T, S, U, R> localized) {
		// LocalRing<T, S> localized = new AbstractLocalizedRing<>(fieldOfFractions,
		// getRing().getIdeal(Collections.singletonList(primeIdealGenerator)));
		PolynomialRing<S> localizedPolynomialRing = AbstractPolynomialRing
				.getPolynomialRing(localized.getLocalizedRing(), numberOfVariables(), getComparator());
		Ideal<Polynomial<S>> localizedIdeal = localizedPolynomialRing.getEmbedding(t, localized.getEmbedding());
		PrimaryDecompositionResult<Polynomial<S>, PolynomialIdeal<S>> decomposition = localizedPolynomialRing
				.primaryDecomposition(localizedIdeal);
		List<PolynomialIdeal<T>> primaries = new ArrayList<>();
		for (PolynomialIdeal<S> primary : decomposition.getPrimaryIdeals()) {
			primaries.add(getEmbedding(primary, localized.getNumerator()));
		}
		List<PolynomialIdeal<T>> radicals = new ArrayList<>();
		for (PolynomialIdeal<S> radical : decomposition.getRadicals()) {
			radicals.add(getEmbedding(radical, localized.getNumerator()));
		}
		return new PrimaryDecompositionResult<>(primaries, radicals);
	}

	private PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionPrincipalIdealRing(
			PolynomialIdeal<T> t, T primeIdealGenerator) {
		if (t.dimension() == 0) {
			if (primeIdealGenerator.equals(getRing().zero())) {
				return primaryDecompositionZeroDimensionalFieldOfFractions(t, getRing().fieldOfFractions());
			}
			return primaryDecompositionZeroDimensionalLocalized(t, primeIdealGenerator,
					getRing().localizeAtIdeal(getRing().getIdeal(Collections.singletonList(primeIdealGenerator))));
		}
		Set<Integer> boundVariables = t.divideOut().boundVariables();
		int variable = 1;
		while (boundVariables.contains(variable)) {
			variable++;
		}
		int[] map = new int[numberOfVariables()];
		int[] inverseMap = new int[numberOfVariables()];
		for (int i = 0; i < numberOfVariables(); i++) {
			if (i < variable - 1) {
				map[i] = i + 1;
				inverseMap[i + 1] = i;
			} else if (i == variable - 1) {
				map[i] = 0;
				inverseMap[0] = i;
			} else {
				map[i] = i;
				inverseMap[i] = i;
			}
		}
		UnivariatePolynomialRing<T> univariate = getRing().getUnivariatePolynomialRing();
		PolynomialRing<Polynomial<T>> newPolynomialRing = AbstractPolynomialRing.getPolynomialRing(univariate,
				numberOfVariables() - 1, getComparator());
		List<Polynomial<Polynomial<T>>> newGenerators = new ArrayList<>();
		for (Polynomial<T> generator : t.generators()) {
			newGenerators.add(unflattenPolynomial(getEmbedding(generator, map), newPolynomialRing, univariate));
		}
		Ideal<Polynomial<T>> newIntersection = newPolynomialRing.getIdeal(newGenerators).intersectToRing();
		PolynomialIdeal<T> saturatedIdeal = t;
		Polynomial<T> newPrimeIdealGenerator = univariate.getEmbedding(primeIdealGenerator);
		for (Polynomial<T> intersectionGenerator : newIntersection.generators()) {
			if (!univariate.isDivisible(intersectionGenerator, newPrimeIdealGenerator)) {
				saturatedIdeal = saturatedIdeal
						.saturate(getEmbedding(intersectionGenerator, new int[] { variable - 1 }));
			}
		}
		if (saturatedIdeal.contains(one())) {
			return new PrimaryDecompositionResult<>(Collections.emptyList(), Collections.emptyList());
		}
		newGenerators.clear();
		for (Polynomial<T> generator : saturatedIdeal.generators()) {
			newGenerators.add(unflattenPolynomial(getEmbedding(generator, map), newPolynomialRing, univariate));
		}
		PolynomialIdeal<Polynomial<T>> newIdeal = newPolynomialRing.getIdeal(newGenerators);
		if (!newIdeal.intersectToRing()
				.equals(univariate.getIdeal(Collections.singletonList(newPrimeIdealGenerator)))) {
			throw new ArithmeticException("Saturation did not work!");
		}
		PrimaryDecompositionResult<Polynomial<Polynomial<T>>, PolynomialIdeal<Polynomial<T>>> decomposition = ((AbstractPolynomialRing<Polynomial<T>>) newPolynomialRing)
				.primaryDecompositionPrincipalIdealRing(newIdeal, newPrimeIdealGenerator);
		List<PolynomialIdeal<T>> primaries = new ArrayList<>();
		List<PolynomialIdeal<T>> radicals = new ArrayList<>();
		for (PolynomialIdeal<Polynomial<T>> primary : decomposition.getPrimaryIdeals()) {
			List<Polynomial<T>> generators = new ArrayList<>();
			for (Polynomial<Polynomial<T>> generator : primary.generators()) {
				generators.add(getEmbedding(flattenPolynomial(generator), inverseMap));
			}
			primaries.add(getIdeal(generators));
		}
		for (PolynomialIdeal<Polynomial<T>> radical : decomposition.getRadicals()) {
			List<Polynomial<T>> generators = new ArrayList<>();
			for (Polynomial<Polynomial<T>> generator : radical.generators()) {
				generators.add(getEmbedding(flattenPolynomial(generator), inverseMap));
			}
			radicals.add(getIdeal(generators));
		}
		Polynomial<T> newDenominator = getEmbedding(((AbstractPolynomialRing<Polynomial<T>>) newPolynomialRing)
				.findDenominator(newIdeal, newPrimeIdealGenerator), new int[] { variable - 1 });
		PolynomialIdeal<T> combined = add(t, getIdeal(Collections.singletonList(newDenominator)));
		if (!combined.contains(one())) {
			PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> newDecomposition = primaryDecompositionPrincipalIdealRing(
					combined, primeIdealGenerator);
			primaries.addAll(newDecomposition.getPrimaryIdeals());
			radicals.addAll(newDecomposition.getRadicals());
		}
		return new PrimaryDecompositionResult<>(primaries, radicals);
	}

	private PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecompositionPrincipalIdealRing(
			PolynomialIdeal<T> t) {
		T denominator = findDenominator(t, getRing().zero());
		PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> decomposition = primaryDecompositionPrincipalIdealRing(
				t, getRing().zero());
		Ideal<T> newDenominatorIdeal = add(t, getIdeal(Collections.singletonList(getEmbedding(denominator))))
				.intersectToRing();
		if (newDenominatorIdeal.generators().size() != 1) {
			throw new ArithmeticException("Not a principal ideal domain");
		}
		T newDenominator = newDenominatorIdeal.generators().get(0);
		if (getRing().isUnit(newDenominator)) {
			return decomposition;
		}
		List<PolynomialIdeal<T>> primaries = new ArrayList<>();
		List<PolynomialIdeal<T>> radicals = new ArrayList<>();
		primaries.addAll(decomposition.getPrimaryIdeals());
		radicals.addAll(decomposition.getRadicals());
		FactorizationResult<T, T> factorization = getRing().uniqueFactorization(newDenominator);
		for (T factor : factorization.primeFactors()) {
			T power = getRing().power(factor, factorization.multiplicity(factor));
			PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> localDecomposition = primaryDecompositionPrincipalIdealRing(
					add(t, getIdeal(Collections.singletonList(getEmbedding(power)))), factor);
			primaries.addAll(localDecomposition.getPrimaryIdeals());
			radicals.addAll(localDecomposition.getRadicals());
		}
		return new PrimaryDecompositionResult<>(primaries, radicals);
	}

	// file:///Users/sophie/Downloads/1-s2.0-S0747717188800403-main.pdf
	@Override
	public PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecomposition(Ideal<Polynomial<T>> t) {
		PolynomialIdeal<T> ideal = (PolynomialIdeal<T>) t;
		if (!getRing().isIntegral() || !getRing().isPrincipalIdealDomain()) {
			throw new UnsupportedOperationException("Not implemented");
		}
		if (ideal.dimension() == 0) {
			Ideal<T> intersectionIdeal = ideal.intersectToRing();
			if (intersectionIdeal.isMaximal()) {
				return primaryDecompositionZeroDimensional(ideal, intersectionIdeal);
			}
		}
		return primaryDecompositionPrincipalIdealRing(ideal);
	}

	@Override
	public ModuloMaximalIdealResult<Polynomial<T>, ?, PolynomialRing<T>, PolynomialIdeal<T>, ?> moduloMaximalIdeal(
			Ideal<Polynomial<T>> ideal) {
		PolynomialIdeal<T> t = (PolynomialIdeal<T>) ideal;
		Ideal<T> intersectToRing = t.intersectToRing();
		return moduloMaximalIdeal(t, getRing().moduloMaximalIdeal(intersectToRing));
	}

	private <S extends Element<S>, R extends Ring<T>, I extends Ideal<T>, F extends Field<S>> ModuloMaximalIdealResult<Polynomial<T>, ?, PolynomialRing<T>, PolynomialIdeal<T>, ?> moduloMaximalIdeal(
			PolynomialIdeal<T> ideal, ModuloMaximalIdealResult<T, S, R, I, F> ringByMaximal) {
		PolynomialRing<S> polynomialRing = AbstractPolynomialRing.getPolynomialRing(ringByMaximal.getField(),
				numberOfVariables(), Monomial.LEX);
		List<Polynomial<S>> reducedGenerators = new ArrayList<>();
		for (Polynomial<T> generator : ideal.generators()) {
			reducedGenerators.add(polynomialRing.getEmbedding(generator, ringByMaximal.getReduction()));
		}
		PolynomialIdeal<S> reducedIdeal = polynomialRing.getIdeal(reducedGenerators);
		return moduloMaximalIdeal(ideal, reducedIdeal, ringByMaximal, polynomialRing,
				ringByMaximal.getField().getExtension(ringByMaximal.getField().getUnivariatePolynomialRing().getVar()));
	}

	private <S extends Element<S>, B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, Ext extends FieldExtension<B, E, Ext>, R extends Ring<T>, I extends Ideal<T>, F extends Field<S>> ModuloMaximalIdealResult<Polynomial<T>, E, PolynomialRing<T>, PolynomialIdeal<T>, Ext> moduloMaximalIdeal(
			PolynomialIdeal<T> ideal, PolynomialIdeal<S> reducedIdeal,
			ModuloMaximalIdealResult<T, S, R, I, F> ringByMaximal, PolynomialRing<S> ringByMaximalPolynomialRing,
			Extension<S, B, E, Ext> fieldAsExtension) {
		PolynomialRing<E> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fieldAsExtension.extension(),
				numberOfVariables(), Monomial.LEX);
		List<Polynomial<E>> reducedGenerators = new ArrayList<>();
		for (Polynomial<S> generator : reducedIdeal.generators()) {
			reducedGenerators.add(polynomialRing.getEmbedding(generator, fieldAsExtension.embeddingMap()));
		}
		PolynomialIdeal<E> extensionIdeal = polynomialRing.getIdeal(reducedGenerators);
		int n = numberOfVariables();
		Ext extension = fieldAsExtension.extension();
		final PolynomialRing<E> originalPolynomialRing = polynomialRing;
		MathMap<Polynomial<T>, Polynomial<E>> reduce = new MathMap<>() {
			@Override
			public Polynomial<E> evaluate(Polynomial<T> t) {
				return originalPolynomialRing.getEmbedding(
						ringByMaximalPolynomialRing.getEmbedding(t, ringByMaximal.getReduction()),
						fieldAsExtension.embeddingMap());
			}
		};
		MathMap<Polynomial<E>, Polynomial<T>> lift = new MathMap<>() {
			@Override
			public Polynomial<T> evaluate(Polynomial<E> t) {
				return getEmbedding(ringByMaximalPolynomialRing.getEmbedding(t, fieldAsExtension.retractionMap()),
						ringByMaximal.getLift());
			}
		};
		while (n > 0) {
			UnivariatePolynomialRing<E> univariateRing = extension.getUnivariatePolynomialRing();
			List<Polynomial<E>> univariateGenerators = new ArrayList<>();
			List<Polynomial<E>> otherGenerators = new ArrayList<>();
			generatorLoop: for (Polynomial<E> generator : extensionIdeal.generators()) {
				for (int i = 1; i < n; i++) {
					if (generator.degree(i) != 0) {
						otherGenerators.add(generator);
						continue generatorLoop;
					}
				}
				univariateGenerators.add(univariateRing.getEmbedding(generator, new int[n]));
			}
			PolynomialIdeal<E> univariateIdeal = univariateRing.getIdeal(univariateGenerators);
			if (univariateIdeal.generators().size() != 1) {
				throw new ArithmeticException("Buchberger failed!");
			}
			FieldEmbedding<B, E, Ext> nextExtension = extension
					.getEmbeddedExtension(univariateRing.toUnivariate(univariateIdeal.generators().get(0)));
			PolynomialRing<E> nextPolynomialRing = AbstractPolynomialRing.getPolynomialRing(nextExtension.getField(),
					n - 1, Monomial.LEX);
			final int nCopy = n;
			MathMap<Polynomial<E>, Polynomial<E>> embedding = new MathMap<>() {
				@Override
				public Polynomial<E> evaluate(Polynomial<E> t) {
					Map<Monomial, E> result = new TreeMap<>();
					for (Monomial m : t.monomials()) {
						int[] exp = Arrays.copyOf(m.exponents(), nCopy - 1);
						E coeff = nextExtension.fromPolynomial(univariateRing.multiply(t.coefficient(m),
								univariateRing.getVarPower(m.exponents()[nCopy - 1])));
						result.put(nextPolynomialRing.getMonomial(exp), coeff);
					}
					return nextPolynomialRing.getPolynomial(result);
				}
			};
			reduce = new ConCatMap<>(reduce, embedding);
			final PolynomialRing<E> oldPolynomialRing = polynomialRing;
			lift = new ConCatMap<Polynomial<E>, Polynomial<E>, Polynomial<T>>(new MathMap<>() {
				@Override
				public Polynomial<E> evaluate(Polynomial<E> t) {
					Map<Monomial, E> result = new TreeMap<>();
					for (Monomial m : t.monomials()) {
						int[] exp = Arrays.copyOf(m.exponents(), nCopy);
						UnivariatePolynomial<E> asPolynomial = nextExtension.asPolynomial(t.coefficient(m));
						for (int i = 0; i <= asPolynomial.degree(); i++) {
							E coeff = asPolynomial.univariateCoefficient(i);
							exp[nCopy - 1] = i;
							Monomial nextM = oldPolynomialRing.getMonomial(exp);
							result.put(nextM, coeff);
						}
					}
					return oldPolynomialRing.getPolynomial(result);
				}
			}, lift);
			List<Polynomial<E>> nextGenerators = new ArrayList<>();
			for (Polynomial<E> generator : otherGenerators) {
				nextGenerators.add(embedding.evaluate(generator));
			}
			PolynomialIdeal<E> nextExtensionIdeal = nextPolynomialRing.getIdeal(nextGenerators);
			polynomialRing = nextPolynomialRing;
			extensionIdeal = nextExtensionIdeal;
			extension = nextExtension.getField();
			n--;
		}
		final MathMap<Polynomial<T>, Polynomial<E>> reduceCopy = reduce;
		final MathMap<Polynomial<E>, Polynomial<T>> liftCopy = lift;
		final PolynomialRing<E> polynomialRingCopy = polynomialRing;
		return new ModuloMaximalIdealResult<>(this, ideal, extension, new MathMap<>() {
			@Override
			public E evaluate(Polynomial<T> t) {
				Polynomial<E> result = reduceCopy.evaluate(t);
				if (result.degree() > 0) {
					throw new ArithmeticException("Not all variables eliminated!");
				}
				return result.leadingCoefficient();
			}
		}, new MathMap<>() {
			@Override
			public Polynomial<T> evaluate(E t) {
				return liftCopy.evaluate(polynomialRingCopy.getEmbedding(t));
			}
		});
	}
//
//	@Override
//	public ModuloMaximalIdealResult<Polynomial<T>, ?> localizeAndModOutMaximalIdeal(Ideal<Polynomial<T>> ideal) {
//		PolynomialIdeal<T> t = (PolynomialIdeal<T>) ideal;
//		Ideal<T> intersected = t.intersectToRing();
//		return localizeAndModOutMaximalIdeal(t, getRing().localizeAndModOutMaximalIdeal(intersected));
//	}
//
//	private <S extends Element<S>> ModuloMaximalIdealResult<Polynomial<T>, ?> localizeAndModOutMaximalIdeal(
//			PolynomialIdeal<T> ideal, ModuloMaximalIdealResult<T, S> ringLocalizedModIdeal) {
//		if (ideal.dimension() == 0) {
//			PolynomialRing<S> localModRing = AbstractPolynomialRing.getPolynomialRing(ringLocalizedModIdeal.getField(),
//					numberOfVariables(), getComparator());
//			PolynomialIdeal<S> t = localModRing.getEmbedding(ideal, ringLocalizedModIdeal.getReduction());
//			return localizeAndModOutMaximalIdealZeroDimensional(ideal, ringLocalizedModIdeal,
//					localModRing.moduloMaximalIdeal(t), localModRing);
//		}
//		int[] map = new int[numberOfVariables()];
//		int[] inverseMap = new int[numberOfVariables()];
//		Set<Integer> boundVariables = ideal.divideOut().boundVariables();
//		int freeIndex = 0;
//		int boundIndex = ideal.dimension();
//		for (int i = 0; i < numberOfVariables(); i++) {
//			if (boundVariables.contains(i + 1)) {
//				map[i] = boundIndex;
//				inverseMap[boundIndex] = i;
//				boundIndex++;
//			} else {
//				map[i] = freeIndex;
//				inverseMap[freeIndex] = freeIndex;
//				freeIndex++;
//			}
//		}
//		TranscendentalFieldExtension<S> transcendentalExtension = new TranscendentalFieldExtension<>(
//				ringLocalizedModIdeal.getField(), ideal.dimension());
//		PolynomialRing<TExt<S>> boundRing = AbstractPolynomialRing.getPolynomialRing(transcendentalExtension,
//				boundVariables.size(), getComparator());
//		PolynomialRing<S> freeRing = transcendentalExtension.polynomialRing();
//		PolynomialRing<S> combinedRing = AbstractPolynomialRing.getPolynomialRing(ringLocalizedModIdeal.getField(),
//				numberOfVariables(), getComparator());
//		PolynomialRing<Polynomial<S>> unflattenRing = AbstractPolynomialRing.getPolynomialRing(freeRing,
//				boundVariables.size(), getComparator());
//		MathMap<Polynomial<T>, Polynomial<TExt<S>>> embedding = new MathMap<>() {
//			@Override
//			public Polynomial<TExt<S>> evaluate(Polynomial<T> t) {
//				Polynomial<S> reordered = combinedRing.getEmbedding(getEmbedding(t, map),
//						ringLocalizedModIdeal.getReduction());
//				Polynomial<Polynomial<S>> unflattened = combinedRing.unflattenPolynomial(reordered, unflattenRing,
//						freeRing);
//				return boundRing.getEmbedding(unflattened, new MathMap<>() {
//					@Override
//					public TExt<S> evaluate(Polynomial<S> t) {
//						return transcendentalExtension.getEmbedding(t);
//					}
//				});
//			}
//		};
//		MathMap<Polynomial<TExt<S>>, Polynomial<T>> asInteger = new MathMap<>() {
//			@Override
//			public Polynomial<T> evaluate(Polynomial<TExt<S>> t) {
//				Polynomial<Polynomial<S>> unflattened = unflattenRing.getEmbedding(t, new MathMap<>() {
//					@Override
//					public Polynomial<S> evaluate(TExt<S> t) {
//						return t.asInteger();
//					}
//				});
//				Polynomial<S> flattened = combinedRing.flattenPolynomial(unflattened);
//				return getEmbedding(getEmbedding(flattened, ringLocalizedModIdeal.getLift()), inverseMap);
//			}
//		};
//		Ideal<Polynomial<TExt<S>>> transcendentalIdeal = boundRing.getIdealEmbedding(ideal, embedding);
//		return localizeAndModOutMaximalIdeal(ideal, embedding, asInteger,
//				boundRing.moduloMaximalIdeal(transcendentalIdeal));
//	}
//
//	private <S extends Element<S>, U extends Element<U>> ModuloMaximalIdealResult<Polynomial<T>, U> localizeAndModOutMaximalIdealZeroDimensional(
//			PolynomialIdeal<T> ideal, ModuloMaximalIdealResult<T, S> ringLocalizedModIdeal,
//			ModuloMaximalIdealResult<Polynomial<S>, U> polynomialRingModIdeal, PolynomialRing<S> polynomialRing) {
//		return new ModuloMaximalIdealResult<>(this, ideal, polynomialRingModIdeal.getField(), new MathMap<>() {
//			@Override
//			public U evaluate(Polynomial<T> t) {
//				return polynomialRingModIdeal.getReduction()
//						.evaluate(polynomialRing.getEmbedding(t, ringLocalizedModIdeal.getReduction()));
//			}
//		}, new MathMap<>() {
//			@Override
//			public Polynomial<T> evaluate(U t) {
//				return getEmbedding(polynomialRingModIdeal.getLift().evaluate(t), ringLocalizedModIdeal.getLift());
//			}
//		});
//	}
//
//	private <S extends Element<S>, U extends Element<U>> ModuloMaximalIdealResult<Polynomial<T>, U> localizeAndModOutMaximalIdeal(
//			PolynomialIdeal<T> ideal, MathMap<Polynomial<T>, Polynomial<TExt<S>>> embedding,
//			MathMap<Polynomial<TExt<S>>, Polynomial<T>> asInteger,
//			ModuloMaximalIdealResult<Polynomial<TExt<S>>, U> transcendentalMod) {
//		return new ModuloMaximalIdealResult<>(this, ideal, transcendentalMod.getField(),
//				new ConCatMap<>(embedding, transcendentalMod.getReduction()),
//				new ConCatMap<>(transcendentalMod.getLift(), asInteger));
//	}

	@Override
	public ModuloIdealResult<Polynomial<T>, ?> moduloIdeal(Ideal<Polynomial<T>> ideal) {
		PolynomialIdeal<T> t = (PolynomialIdeal<T>) ideal;
		return moduloIdeal(t, getRing().moduloIdeal(t.intersectToRing()));
	}

	private <S extends Element<S>> ModuloIdealResult<Polynomial<T>, CoordinateRingElement<S>> moduloIdeal(
			PolynomialIdeal<T> ideal, ModuloIdealResult<T, S> moduloRing) {
		PolynomialRing<S> moduloPolynomialRing = AbstractPolynomialRing.getPolynomialRing(moduloRing.getQuotientRing(),
				numberOfVariables(), getComparator());
		CoordinateRing<S> result = moduloPolynomialRing.getEmbedding(ideal, moduloRing.getReduction()).divideOut();
		return new ModuloIdealResult<>(this, ideal, result, new MathMap<>() {

			@Override
			public CoordinateRingElement<S> evaluate(Polynomial<T> t) {
				return result.getEmbedding(moduloPolynomialRing.getEmbedding(t, moduloRing.getReduction()));
			}
		}, new MathMap<>() {

			@Override
			public Polynomial<T> evaluate(CoordinateRingElement<S> t) {
				return getEmbedding(t.getElement(), moduloRing.getLift());
			}
		});
	}

	@Override
	public LocalizeResult<Polynomial<T>, ?, ?, ?> localizeAtIdeal(Ideal<Polynomial<T>> primeIdeal) {
		PolynomialIdeal<T> ideal = (PolynomialIdeal<T>) primeIdeal;
		return localizeAtIdeal(ideal, getRing().localizeAtIdeal(ideal.intersectToRing()));
	}

	private <S extends Element<S>, U extends Element<U>, R extends Element<R>> LocalizeResult<Polynomial<T>, LocalizedElement<S>, LocalizedElement<S>, S> localizeAtIdeal(
			PolynomialIdeal<T> primeIdeal, LocalizeResult<T, S, U, R> ringLocalization) {
		PolynomialRing<S> polynomialRing = AbstractPolynomialRing.getPolynomialRing(ringLocalization.getLocalizedRing(),
				numberOfVariables(), getComparator());
		List<Polynomial<S>> generators = new ArrayList<>();
		for (Polynomial<T> generator : primeIdeal.generators()) {
			if (generator.degree() > 0) {
				generators.add(polynomialRing.getEmbedding(generator, ringLocalization.getEmbedding()));
			}
		}
		PolynomialIdeal<S> localizedIdeal = polynomialRing.getIdeal(generators);
		CoordinateRing<S> coordinateRing = polynomialRing.getZeroIdeal().divideOut();
		// TODO: Fix
		LocalizedCoordinateRing<S> result = new LocalizedCoordinateRing<>(
				(Field<S>) ringLocalization.getLocalizedRing(), coordinateRing,
				coordinateRing.getIdeal(localizedIdeal));
		MathMap<LocalizedElement<S>, S> denominator = new MathMap<>() {

			@Override
			public S evaluate(LocalizedElement<S> t) {
				Polynomial<S> numerator = t.asPolynomialFraction().getNumerator();
				T denominator = getRing().one();
				for (Monomial m : numerator.monomials()) {
					denominator = getRing().lcm(denominator,
							ringLocalization.getDenominator().evaluate(numerator.coefficient(m)));
				}
				Polynomial<S> denom = t.asPolynomialFraction().getDenominator();
				for (Monomial m : denom.monomials()) {
					denominator = getRing().lcm(denominator,
							ringLocalization.getDenominator().evaluate(denom.coefficient(m)));
				}
				return ringLocalization.getEmbedding().evaluate(denominator);
			}
		};
		return new LocalizeResult<>(this, primeIdeal, result.ringOfIntegers(), new MathMap<>() {
			@Override
			public LocalizedElement<S> evaluate(Polynomial<T> t) {
				return result.getEmbedding(polynomialRing.getEmbedding(t, ringLocalization.getEmbedding()));
			}
		}, new MathMap<>() {
			@Override
			public Polynomial<T> evaluate(LocalizedElement<S> t) {
				Polynomial<S> numerator = t.asPolynomialFraction().getNumerator();
				numerator = polynomialRing.multiply(denominator.evaluate(t), numerator);
				return getEmbedding(numerator, ringLocalization.getAsInteger());
			}
		}, new MathMap<>() {
			@Override
			public Polynomial<T> evaluate(LocalizedElement<S> t) {
				Polynomial<S> denom = t.asPolynomialFraction().getNumerator();
				denom = polynomialRing.multiply(denominator.evaluate(t), denom);
				return getEmbedding(denom, ringLocalization.getAsInteger());
			}
		}, new MathMap<>() {
			@Override
			public Polynomial<T> evaluate(LocalizedElement<S> t) {
				return getEmbedding(t.asPolynomialFraction().asInteger(), ringLocalization.getAsInteger());
			}
		});
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
	public GroebnerBasis<T> buchberger(List<Polynomial<T>> generators, boolean computeExpressionsAndSyzygies) {
		List<Polynomial<T>> groebner = new ArrayList<>();
		List<List<Polynomial<T>>> expressions = new ArrayList<>();
		List<Vector<Polynomial<T>>> syzygies = new ArrayList<>();
		groebner.addAll(generators);
		for (int i = 0; i < generators.size(); i++) {
			if (!generators.get(i).getPolynomialRing().equals(this)) {
				throw new ArithmeticException("Wrong polynomial ring!");
			}
			List<Polynomial<T>> generatorExpression = new ArrayList<>();
			for (int j = 0; j < generators.size(); j++) {
				generatorExpression.add(zero());
			}
			generatorExpression.set(i, one());
			expressions.add(generatorExpression);
		}
		GroebnerBasis<T> basis = reduceBasis(new GroebnerBasis<>(groebner, expressions, syzygies),
				computeExpressionsAndSyzygies);
		groebner = basis.getBasis();
		expressions = basis.getExpression();
		syzygies = basis.getSyzygies();
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
					if (computeExpressionsAndSyzygies) {
						for (int k = 0; k < generators.size(); k++) {
							newGeneratorExpression.add(subtract(multiply(pli, expressions.get(i).get(k)),
									this.multiply(plj, expressions.get(j).get(k))));
						}
					}
					ReduceAndExpressResult<T> re = reduceAndExpress(newGenerator, newGeneratorExpression, groebner,
							expressions, computeExpressionsAndSyzygies);
					newGenerator = re.getReduced();
					newGeneratorExpression = re.getExpression();
					re = reduceAndExpress(newGenerator, newGeneratorExpression, newGroebner, newExpressions,
							computeExpressionsAndSyzygies);
					newGenerator = re.getReduced();
					newGeneratorExpression = re.getExpression();
					if (!newGenerator.equals(this.zero())) {
						newGroebner.add(newGenerator);
						newExpressions.add(newGeneratorExpression);
					} else if (computeExpressionsAndSyzygies) {
						Vector<Polynomial<T>> syzygy = new Vector<>(newGeneratorExpression);
						syzygies.add(syzygy);
					}
					if (!getRing().isEuclidean() || getRing().krullDimension() == 0) {
						continue;
					}
					ExtendedEuclideanResult<T> extendedEuclid = getRing().extendedEuclidean(ci, cj);
					if (getRing().isDivisible(extendedEuclid.getGcd(), ci)
							|| getRing().isDivisible(extendedEuclid.getGcd(), cj)) {
						continue;
					}
					newGenerator = add(multiply(extendedEuclid.getCoeff1(), li, pi),
							multiply(extendedEuclid.getCoeff2(), lj, pj));
					newGeneratorExpression = new ArrayList<>();
					if (computeExpressionsAndSyzygies) {
						for (int k = 0; k < generators.size(); k++) {
							newGeneratorExpression
									.add(add(multiply(extendedEuclid.getCoeff1(), li, expressions.get(i).get(k)),
											this.multiply(extendedEuclid.getCoeff2(), lj, expressions.get(j).get(k))));

						}
					}
					re = reduceAndExpress(newGenerator, newGeneratorExpression, groebner, expressions,
							computeExpressionsAndSyzygies);
					newGenerator = re.getReduced();
					newGeneratorExpression = re.getExpression();
					re = reduceAndExpress(newGenerator, newGeneratorExpression, newGroebner, newExpressions,
							computeExpressionsAndSyzygies);
					newGenerator = re.getReduced();
					newGeneratorExpression = re.getExpression();
					if (!newGenerator.equals(this.zero())) {
						newGroebner.add(newGenerator);
						newExpressions.add(newGeneratorExpression);
					} else if (computeExpressionsAndSyzygies) {
						Vector<Polynomial<T>> syzygy = new Vector<>(newGeneratorExpression);
						syzygies.add(syzygy);
					}
				}
			}
			newGroebner.addAll(groebner);
			newExpressions.addAll(expressions);
			GroebnerBasis<T> b = reduceBasis(new GroebnerBasis<>(newGroebner, newExpressions, Collections.emptyList()),
					computeExpressionsAndSyzygies);
			newGroebner = b.getBasis();
			newExpressions = b.getExpression();
			if (computeExpressionsAndSyzygies) {
				syzygies.addAll(b.getSyzygies());
				TreeSet<Vector<Polynomial<T>>> syzygiesSet = new TreeSet<>();
				syzygiesSet.addAll(syzygies);
				syzygies.clear();
				syzygies.addAll(syzygiesSet);
			}
			if (groebner.equals(newGroebner)) {
				break;
			}
			groebner = newGroebner;
			expressions = newExpressions;
		}
		return new GroebnerBasis<>(groebner, expressions, syzygies);
	}

	@Override
	public GroebnerBasis<T> reduceBasis(GroebnerBasis<T> basis, boolean computeExpressionsAndSyzygies) {
		if (basis.getBasis().isEmpty()) {
			return basis;
		}
		List<Vector<Polynomial<T>>> syzygies = new ArrayList<>();
		if (computeExpressionsAndSyzygies) {
			syzygies.addAll(basis.getSyzygies());
		}
		List<ReduceAndExpressResult<T>> list = new ArrayList<>();
		for (int i = 0; i < basis.getBasis().size(); i++) {
			if (basis.getBasis().get(i).equals(zero())) {
				if (computeExpressionsAndSyzygies) {
					Vector<Polynomial<T>> syzygy = new Vector<>(basis.getExpression().get(i));
					syzygies.add(syzygy);
				}
			} else {
				list.add(new ReduceAndExpressResult<>(basis.getBasis().get(i), basis.getExpression().get(i), false));
			}
		}
		Collections.sort(list, Comparator.reverseOrder());
		boolean modified = true;
		while (modified) {
			modified = false;
			candidateLoop: for (int i = 0; i < list.size(); i++) {
				ReduceAndExpressResult<T> polynomial = list.remove(i);
				ReduceAndExpressResult<T> reduce = this.reduceAndExpress(polynomial.getReduced(),
						polynomial.getExpression(), list, computeExpressionsAndSyzygies);
				modified = reduce.isModified();
				if (modified) {
					if (reduce.getReduced().equals(zero())) {
						if (computeExpressionsAndSyzygies) {
							Vector<Polynomial<T>> syzygy = new Vector<>(reduce.getExpression());
							syzygies.add(syzygy);
						}
					} else {
						for (int j = 0; j < list.size(); j++) {
							if (reduce.compareTo(list.get(j)) > 0) {
								list.add(j, reduce);
								break candidateLoop;
							}
						}
						list.add(reduce);
					}
					break;
				} else {
					list.add(i, reduce);
				}
			}
		}
		List<Polynomial<T>> reducedBasis = new ArrayList<>();
		List<List<Polynomial<T>>> reducedExpressions = new ArrayList<>();
		for (int i = 0; i < list.size(); i++) {
			reducedBasis.add(list.get(i).getReduced());
			reducedExpressions.add(list.get(i).getExpression());
		}
		return new GroebnerBasis<>(reducedBasis, reducedExpressions, syzygies);
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
		return getRing() + "[" + String.join(",", variableNames) + "]";
	}

	public Monomial parseMonomial(PeekableReader in) throws IOException {
		if (variableTrie == null) {
			variableTrie = new Trie(variableNames);
			if (!variableTrie.isPrefixFree()) {
				throw new IOException("Set of variables not prefix free!");
			}
			variableMap = new TreeMap<>();
			for (int i = 0; i < variableNames.length; i++) {
				variableMap.put(variableNames[i], i);
			}
		}
		Integers z = Integers.z();
		int variablesParsed = 0;
		int[] exp = new int[numberOfVariables()];
		Node variable = variableTrie.getRootNode();
		int length = 0;
		while (true) {
			int result = in.peek(length);
			length++;
			if (result < 0) {
				break;
			}
			Node nextVariable = variable.child((char) result);
			if (nextVariable == null) {
				break;
			}
			if (nextVariable.isValid()) {
				String var = nextVariable.reconstruct();
				in.skip(length);
				length = 0;
				variable = variableTrie.getRootNode();
				int exponent = 1;
				if (in.peek() == '^') {
					in.skip(1);
					exponent = z.parse(in).intValueExact();
				}
				exp[variableMap.get(var)] = exponent;
				variablesParsed++;
				if (in.peek() == '*') {
					in.skip(1);
					continue;
				}
				break;
			}
			variable = nextVariable;
		}
		if (variablesParsed == 0) {
			throw new IOException("No monomial found!");
		}
		return getMonomial(exp);
	}

	@Override
	public Polynomial<T> parse(PeekableReader reader) throws IOException {
		int result = reader.peek();
		boolean brackets = false;
		if (result == '(') {
			brackets = true;
			reader.skip(1);
		}
		Ring<T> r = getRing();
		boolean plusFound = true;
		Map<Monomial, T> parsed = new TreeMap<>();
		while (true) {
			int ch = reader.peek();
			if (ch < 0) {
				if (plusFound) {
					throw new IOException("Malformed polynomial!");
				}
				break;
			}
			if (ch == ' ') {
				reader.skip(1);
				continue;
			}
			if (!plusFound) {
				if (ch == '+') {
					plusFound = true;
					reader.skip(1);
					continue;
				}
				break;
			}
			T coeff;
			Monomial m;
			try {
				coeff = r.parse(reader);
				if (reader.peek() == '*') {
					reader.skip(1);
				} else {
					parsed.put(getMonomial(new int[numberOfVariables()]), coeff);
					plusFound = false;
					continue;
				}
			} catch (IOException e) {
				coeff = r.one();
			}
			try {
				m = parseMonomial(reader);
			} catch (IOException e) {
				throw new IOException("malformed polynomial!", e);
			}
			parsed.put(m, coeff);
			plusFound = false;
		}
		if (parsed.size() == 0) {
			throw new IOException("malformed polynomial!");
		}
		if (brackets) {
			int ch = reader.peek();
			if (ch < 0 || ch != ')') {
				throw new IOException("mismatched brackets!");
			}
			reader.skip(1);
		}
		return getPolynomial(parsed);
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
				if (evaluated.degree(1) != degreeX || evaluated.degree(2) != degreeY) {
					continue mainLoop;
				}
			}
			List<T> evaluationPoint = new ArrayList<>();
			evaluationPoint.add(null);
			evaluationPoint.add(null);
			evaluationPoint.addAll(pointList);
			Polynomial<T> evaluated = bivariate.getEmbedding(partiallyEvaluate(t, evaluationPoint));
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
			throw new UnsupportedOperationException("Not yet further implemented");
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
		if (numberOfVariables() == 0) {
			T unit = getRing().projectToUnit(t.leadingCoefficient());
			return new FactorizationResult<>(unit, SingletonSortedMap.map(divideScalar(t, unit), 1));
		}
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
			return new FactorizationResult<Polynomial<T>, Polynomial<T>>(getEmbedding(factors.getUnit()),
					factors.factorMap());
		}
		if (!getRing().isUniqueFactorizationDomain()) {
			throw new ArithmeticException("Not a UFD!");
		}
		if (getRing().krullDimension() != 0) {
			FieldOfFractions<T> fieldOfFractions = new FieldOfFractions<>(getRing());
			T content = content(t);
			FactorizationResult<T, T> contentFactors = getRing().uniqueFactorization(content);
			t = contentFree(t);
			PolynomialRing<Fraction<T>> fractionPolynomialRing = getPolynomialRing(fieldOfFractions,
					numberOfVariables(), getComparator());
			Polynomial<Fraction<T>> overField = fractionPolynomialRing.getEmbedding(t, new MathMap<>() {
				@Override
				public Fraction<T> evaluate(T t) {
					return fieldOfFractions.getEmbedding(t);
				}
			});
			FactorizationResult<Polynomial<Fraction<T>>, Polynomial<Fraction<T>>> factors = fractionPolynomialRing
					.uniqueFactorization(overField);
			Fraction<T> unit = fieldOfFractions.multiply(fieldOfFractions.getEmbedding(contentFactors.getUnit()),
					factors.getUnit().leadingCoefficient());
			SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
			for (Polynomial<Fraction<T>> factor : factors.primeFactors()) {
				int multiplicity = factors.multiplicity(factor);
				factor = fractionPolynomialRing.scalarMultiply(fieldOfFractions.getEmbedding(content), factor);
				Polynomial<T> intFactor = getEmbedding(factor, new MathMap<>() {
					@Override
					public T evaluate(Fraction<T> t) {
						return t.asInteger();
					}
				});
				unit = fieldOfFractions.multiply(unit, fieldOfFractions.getEmbedding(content(intFactor)));
				intFactor = contentFree(intFactor);
				result.put(intFactor, multiplicity);
			}
			for (T prime : contentFactors.primeFactors()) {
				result.put(getEmbedding(prime), contentFactors.multiplicity(prime));
			}
			return new FactorizationResult<>(getEmbedding(unit.asInteger()), result);
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
		if (polynomials.isEmpty()) {
			throw new InfinityException();
		}
		PolynomialIdeal<T> ideal = getIdeal(polynomials);
		return solve(ideal);
	}

	@Override
	public List<AffinePoint<T>> solve(PolynomialIdeal<T> ideal) {
		Field<T> field = (Field<T>) getRing();
		if (ideal.contains(one())) {
			return Collections.emptyList();
		}
		if (ideal.dimension() != 0) {
			throw new InfinityException();
		}
		if (numberOfVariables() == 1) {
			Set<T> roots = new TreeSet<>();
			roots.addAll(field.roots(ideal.generators().get(0)).keySet());
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
			for (Polynomial<T> polynomial : ideal.generators()) {
				embeded.add(revlexRing.getEmbedding(polynomial, map));
			}
			return revlexRing.solve(embeded);
		}
		int[] mapShift = new int[numberOfVariables()];
		int[] mapFirst = new int[numberOfVariables()];
		mapShift[0] = -1;
		mapFirst[0] = 0;
		for (int i = 0; i < numberOfVariables() - 1; i++) {
			mapShift[i + 1] = i;
			mapFirst[i + 1] = -1;
		}
		int numberOfPolynomials = ideal.generators().size();
		UnivariatePolynomialRing<T> univariate = field.getUnivariatePolynomialRing();
		PolynomialRing<T> shifted = AbstractPolynomialRing.getPolynomialRing(field, numberOfVariables() - 1,
				Monomial.REVLEX);
		UnivariatePolynomial<T> last = univariate.getEmbedding(ideal.generators().get(numberOfPolynomials - 1),
				mapFirst);
		List<AffinePoint<T>> result = new ArrayList<>();
		for (T root : field.roots(last).keySet()) {
			List<T> eval = new ArrayList<>();
			eval.add(root);
			for (int i = 1; i < numberOfVariables(); i++) {
				eval.add(null);
			}
			List<Polynomial<T>> evaluated = new ArrayList<>();
			for (Polynomial<T> generator : ideal.generators()) {
				evaluated.add(shifted.getEmbedding(partiallyEvaluate(generator, eval), mapShift));
			}
			for (AffinePoint<T> partial : shifted.solve(evaluated)) {
				List<T> point = new ArrayList<>();
				point.add(root);
				point.addAll(partial.getCoords());
				result.add(new AffinePoint<>(field, point));
			}
		}
		return result;
	}

	private Polynomial<T> rowMultiply(Vector<Polynomial<T>> row, Vector<Polynomial<T>> column) {
		if (row.dimension() != column.dimension()) {
			throw new ArithmeticException("Mismatched dimensions");
		}
		Polynomial<T> result = zero();
		for (int i = 0; i < row.dimension(); i++) {
			result = add(multiply(row.get(i + 1), column.get(i + 1)), result);
		}
		return result;
	}

	private List<Polynomial<T>> rowMultiplyMatrix(Vector<Polynomial<T>> row, Matrix<Polynomial<T>> columns) {
		List<Polynomial<T>> result = new ArrayList<>();
		for (int i = 0; i < columns.columns(); i++) {
			result.add(rowMultiply(row, columns.column(i + 1)));
		}
		return result;
	}

	@Override
	public boolean isSubModuleMember(MatrixModule<Polynomial<T>> module, Matrix<Polynomial<T>> m,
			Vector<Polynomial<T>> b) {
		Vector<Polynomial<T>> asSubModuleMember = module.domain().zero();
		Matrix<Polynomial<T>> syzygiesMatrix = module.domainAlgebra().one();
		for (int i = 0; i < m.rows(); i++) {
			Vector<Polynomial<T>> row = m.row(i + 1);
			IdealResult<Polynomial<T>, PolynomialIdeal<T>> ideal = getIdealWithTransforms(
					rowMultiplyMatrix(row, syzygiesMatrix));
			Polynomial<T> adjusted = subtract(b.get(i + 1), rowMultiply(row, asSubModuleMember));
			if (!ideal.getIdeal().contains(adjusted)) {
				return false;
			}
			List<Polynomial<T>> idealResult = ideal.getIdeal().generate(adjusted);
			for (int j = 0; j < ideal.getIdeal().generators().size(); j++) {
				asSubModuleMember = module
						.domain().add(
								module.domain().scalarMultiply(idealResult.get(j),
										module.multiply(syzygiesMatrix,
												new Vector<>(ideal.getGeneratorExpressions().get(j)))),
								asSubModuleMember);
			}
			if (ideal.getSyzygies().size() > 0) {
				Matrix<Polynomial<T>> rowSyzygiesMatrix = Matrix.fromColumns(ideal.getSyzygies());
				syzygiesMatrix = module.multiply(syzygiesMatrix, rowSyzygiesMatrix);
			} else {
				syzygiesMatrix = module.domainAlgebra().zero();
			}
		}
		return true;
	}

	@Override
	public Vector<Polynomial<T>> asSubModuleMember(MatrixModule<Polynomial<T>> module, Matrix<Polynomial<T>> m,
			Vector<Polynomial<T>> b) {
		Vector<Polynomial<T>> asSubModuleMember = module.domain().zero();
		Matrix<Polynomial<T>> syzygiesMatrix = module.domainAlgebra().one();
		for (int i = 0; i < m.rows(); i++) {
			Vector<Polynomial<T>> row = m.row(i + 1);
			IdealResult<Polynomial<T>, PolynomialIdeal<T>> ideal = getIdealWithTransforms(
					rowMultiplyMatrix(row, syzygiesMatrix));
			Polynomial<T> adjusted = subtract(b.get(i + 1), rowMultiply(row, asSubModuleMember));
			List<Polynomial<T>> idealResult = ideal.getIdeal().generate(adjusted);
			for (int j = 0; j < ideal.getIdeal().generators().size(); j++) {
				asSubModuleMember = module
						.domain().add(
								module.domain().scalarMultiply(idealResult.get(j),
										module.multiply(syzygiesMatrix,
												new Vector<>(ideal.getGeneratorExpressions().get(j)))),
								asSubModuleMember);
			}
			if (ideal.getSyzygies().size() > 0) {
				Matrix<Polynomial<T>> rowSyzygiesMatrix = Matrix.fromColumns(ideal.getSyzygies());
				syzygiesMatrix = module.multiply(syzygiesMatrix, rowSyzygiesMatrix);
			} else {
				syzygiesMatrix = module.domainAlgebra().zero();
			}
		}
		return asSubModuleMember;
	}

	@Override
	public List<Vector<Polynomial<T>>> syzygyProblem(MatrixModule<Polynomial<T>> module, Matrix<Polynomial<T>> m) {
		Matrix<Polynomial<T>> syzygiesMatrix = module.domainAlgebra().one();
		for (int i = 0; i < m.rows(); i++) {
			List<Vector<Polynomial<T>>> rowSyzygies = getIdealWithTransforms(
					rowMultiplyMatrix(m.row(i + 1), syzygiesMatrix)).getSyzygies();
			if (rowSyzygies.size() > 0) {
				Matrix<Polynomial<T>> rowSyzygiesMatrix = Matrix.fromColumns(rowSyzygies);
				syzygiesMatrix = module.multiply(syzygiesMatrix, rowSyzygiesMatrix);
			} else {
				return Collections.emptyList();
			}
		}
		return simplifySubModuleGenerators(syzygiesMatrix.getModule(this), syzygiesMatrix);
	}

	@Override
	public List<Vector<Polynomial<T>>> simplifySubModuleGenerators(MatrixModule<Polynomial<T>> module,
			Matrix<Polynomial<T>> m) {
		List<Vector<Polynomial<T>>> columns = new ArrayList<>();
		for (int i = 0; i < m.columns(); i++) {
			columns.add(m.column(i + 1));
		}
		Collections.sort(columns, new Comparator<Vector<Polynomial<T>>>() {
			@Override
			public int compare(Vector<Polynomial<T>> o1, Vector<Polynomial<T>> o2) {
				int maxDegree1 = -1;
				int maxDegree2 = -1;
				for (Polynomial<T> t : o1.asList()) {
					if (maxDegree1 < t.degree()) {
						maxDegree1 = t.degree();
					}
				}
				for (Polynomial<T> t : o2.asList()) {
					if (maxDegree2 < t.degree()) {
						maxDegree2 = t.degree();
					}
				}
				if (maxDegree1 != maxDegree2) {
					return maxDegree1 - maxDegree2;
				}
				return o1.compareTo(o2);
			}
		});
		List<Vector<Polynomial<T>>> result = new ArrayList<>();
		MatrixModule<Polynomial<T>> matrixModule = null;
		for (int i = 0; i < m.columns(); i++) {
			Vector<Polynomial<T>> column = columns.get(i);
			if (!column.equals(module.codomain().zero())) {
				if (result.isEmpty() || !isSubModuleMember(matrixModule, Matrix.fromColumns(result), column)) {
					result.add(column);
					matrixModule = new MatrixModule<>(this, m.rows(), result.size());
				}
			}
		}
		return result;
	}

}
