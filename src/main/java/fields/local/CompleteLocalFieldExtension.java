package fields.local;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.FieldAutomorphism;
import fields.helper.FieldEmbedding;
import fields.helper.GaloisGroup;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;
import fields.interfaces.LocalField;
import fields.interfaces.LocalRing;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.CompleteLocalFieldExtension.Ext;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import util.Identity;
import util.Pair;

public class CompleteLocalFieldExtension<B extends Element<B>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>>
		extends AbstractFieldExtension<B, Ext<B>, CompleteLocalFieldExtension<B, R, RE, RFE>>
		implements LocalField<Ext<B>, RE> {
	private LocalField<B, R> base;
	private LocalRingExtension<B, Ext<B>, CompleteLocalFieldExtension<B, R, RE, RFE>, R, RE, RFE> localRing;
	private FieldEmbedding<B, Ext<B>, CompleteLocalFieldExtension<B, R, RE, RFE>> maximalUnramifiedExtension;
	private RFE reduction;
	private Ext<B> reductionGenerator;
	private Ext<B> uniformizer;
	private int residueDegree;
	private int ramificationIndex;
	private Matrix<B> alphaToReductionGeneratorUniformizerBaseChange;
	private Map<Ext<B>, List<Ext<B>>> conjugates;
	private LocalField<B, R> highAccuracyBase;
	private UnivariatePolynomialRing<B> highAccuracyPolynomials;
	private MatrixAlgebra<B> highAccuracyAlgebra;
	private Reals r = Reals.r(1024);

	public static class Ext<B extends Element<B>> extends AbstractElement<Ext<B>>
			implements AlgebraicExtensionElement<B, Ext<B>> {
		private UnivariatePolynomial<B> e;

		public Ext(UnivariatePolynomial<B> e) {
			this.e = e;
		}

		@Override
		public int compareTo(Ext<B> o) {
			return this.e.compareTo(o.e);
		}

		@Override
		public UnivariatePolynomial<B> asPolynomial() {
			return e;
		}

		@Override
		public String toString() {
			return e.toString("α", true);
		}
	}

	public CompleteLocalFieldExtension(UnivariatePolynomial<B> minimalPolynomial, LocalField<B, R> base,
			RFE reductionExtensionSample) {
		super(minimalPolynomial, base);
		this.base = base;
		this.conjugates = new TreeMap<>();
		this.highAccuracyBase = base.withAccuracy(base.getAccuracy() * degree());
		this.highAccuracyPolynomials = highAccuracyBase.getUnivariatePolynomialRing();
		this.highAccuracyAlgebra = new FreeModule<>(highAccuracyBase, degree()).matrixAlgebra();
		if (base.ringOfIntegers().hasIrreducibleGoodReduction(minimalPolynomial)) {
			initUnramified(reductionExtensionSample);
			return;
		}
		initGeneric(reductionExtensionSample);
	}

	private void initUnramified(RFE reductionExtensionSample) {
		UnivariatePolynomial<R> reducedMipo = base.ringOfIntegers().reduceUnivariatePolynomial(minimalPolynomial());
		this.reduction = reductionExtensionSample.makeExtension(reducedMipo);
		this.maximalUnramifiedExtension = new FieldEmbedding<>(this);
		this.ramificationIndex = 1;
		this.residueDegree = degree();
		this.reductionGenerator = alpha();
		this.uniformizer = getEmbedding(base.uniformizer());
		this.alphaToReductionGeneratorUniformizerBaseChange = matrixAlgebra().one();
		this.localRing = null;//new LocalRingExtension<>(this, base, reduction);
	}

	@SuppressWarnings("unchecked")
	private void initGeneric(RFE reductionExtensionSample) {
		UnivariatePolynomialRing<Ext<B>> polynomials = getUnivariatePolynomialRing();
		LocalRing<B, R> baseRing = base.ringOfIntegers();

		this.reduction = reductionExtensionSample
				.makeExtension(base.reduction().getUnivariatePolynomialRing().getVar());
		this.maximalUnramifiedExtension = new FieldEmbedding<>(this,
				makeExtension(base.getUnivariatePolynomialRing().getVar()), zero());
		this.reductionGenerator = maximalUnramifiedExtension.getEmbeddedAlpha();
		this.localRing = null;//new LocalRingExtension<>(this, base, reduction);
		CompleteLocalFieldExtension<B, R, RE, RFE> unramified = maximalUnramifiedExtension.getEmbeddedField();
		LocalRingExtension<B, Ext<B>, CompleteLocalFieldExtension<B, R, RE, RFE>, R, RE, RFE> unramifiedRing = unramified
				.ringOfIntegers();
		UnivariatePolynomial<Ext<B>> minimalPolynomialOverUnramified = polynomials.getEmbedding(minimalPolynomial(),
				getEmbeddingMap());
		Ext<B> alpha = alpha();

		while (true) {
			if (unramifiedRing.isEisenstein(minimalPolynomialOverUnramified)) {
				this.ramificationIndex = minimalPolynomialOverUnramified.degree();
				this.residueDegree = reduction.degree();
				this.uniformizer = alpha;
				List<Vector<B>> reductionGeneratorUniformizerBase = new ArrayList<>();
				Ext<B> uniformizerPower = one();
				for (int i = 0; i < ramificationIndex; i++) {
					Ext<B> reductionGeneratorPower = one();
					for (int j = 0; j < residueDegree; j++) {
						reductionGeneratorUniformizerBase
								.add(asVector(multiply(uniformizerPower, reductionGeneratorPower)));
						reductionGeneratorPower = multiply(reductionGeneratorPower, reductionGenerator);
					}
					uniformizerPower = multiply(uniformizerPower, uniformizer);
				}
				this.alphaToReductionGeneratorUniformizerBaseChange = matrixAlgebra()
						.inverse(Matrix.fromColumns(reductionGeneratorUniformizerBase));
				return;
			}
			if (minimalPolynomialOverUnramified.degree() == 1) {
				this.ramificationIndex = 1;
				this.residueDegree = reduction.degree();
				this.uniformizer = getEmbedding(base.uniformizer());
				this.alphaToReductionGeneratorUniformizerBaseChange = matrixAlgebra().one();
				return;
			}
			Pair<UnivariatePolynomial<Ext<B>>, UnivariatePolynomial<Ext<B>>> integral = unramifiedRing
					.integralPolynomial(minimalPolynomialOverUnramified,
							unramified.withAccuracy(base.getAccuracy() * degree()).getUnivariatePolynomialRing());
			minimalPolynomialOverUnramified = unramifiedRing.roundUnivariatePolynomial(integral.getFirst(),
					base.getAccuracy());
			alpha = polynomials.evaluate(
					polynomials.getEmbedding(integral.getSecond(), maximalUnramifiedExtension.getEmbeddingMap()),
					alpha);
//			if (!polynomials.evaluate(polynomials.getEmbedding(minimalPolynomialOverUnramified,
//					maximalUnramifiedExtension.getEmbeddingMap()), alpha).equals(zero())) {
//				System.err.println("Alpha not zero of minimalPolynomialOverUnramified! " + polynomials.evaluate(polynomials.getEmbedding(minimalPolynomialOverUnramified,
//						maximalUnramifiedExtension.getEmbeddingMap()), alpha));
//			}
			UnivariatePolynomial<RE> reducedMinimalPolynomialOverUnramified = unramifiedRing
					.reduceUnivariatePolynomial(minimalPolynomialOverUnramified);
			FactorizationResult<Polynomial<RE>, RE> reducedFactors = unramified.reduction()
					.factorization(reducedMinimalPolynomialOverUnramified);
			if (!reducedFactors.isIrreducible()) {
				throw new ArithmeticException("minimal polynomial not irreducible");
			}
			UnivariatePolynomial<RE> reduced = unramified.reduction().getUnivariatePolynomialRing()
					.toUnivariate(reducedFactors.firstPrimeFactor());
			Value constantValuation = unramified.valuation(minimalPolynomialOverUnramified.univariateCoefficient(0));
			if (reduced.degree() != 1) {
				Extension<RE, R, RE, RFE> reductionExtension = reduction.getExtension(reduced);
				this.reduction = reductionExtension.extension();
				this.residueDegree = reduction.degree();
				this.ramificationIndex = degree() / residueDegree;
				this.localRing = null;//new LocalRingExtension<>(this, base, reduction);
				UnivariatePolynomial<B> lifted = baseRing.liftUnivariatePolynomial(reduction.minimalPolynomial());
				unramified = makeExtension(lifted);
				Vector<RE> asReducedVector = reductionExtension.asVectorMap().evaluate(reduction.alpha());
				Ext<B> power = one();
				Ext<B> initial = zero();
				for (int i = 0; i < asReducedVector.dimension(); i++) {
					initial = add(initial, multiply(power,
							maximalUnramifiedExtension.getEmbedding(unramified.lift(asReducedVector.get(i + 1)))));
					power = multiply(power, alpha);
				}
				this.reductionGenerator = localRing.henselLiftWithInitialLift(
						polynomials.getEmbedding(lifted, getEmbeddingMap()), initial,
						2 * base.getAccuracy() * degree());
				this.maximalUnramifiedExtension = new FieldEmbedding<>(this, unramified, reductionGenerator);
				minimalPolynomialOverUnramified = maximalUnramifiedExtension.minimalPolynomial();
				unramifiedRing = unramified.ringOfIntegers();
				alpha = alpha();
			} else if (constantValuation.equals(Value.ZERO)) {
				Ext<B> root = unramifiedRing.lift(unramified.reduction.roots(reduced).keySet().iterator().next());
				Ext<B> embeddedRoot = maximalUnramifiedExtension.getEmbedding(root);
				alpha = subtract(alpha, embeddedRoot);
				UnivariatePolynomialRing<Ext<B>> unramifiedPolynomials = unramified.getUnivariatePolynomialRing();
				minimalPolynomialOverUnramified = unramifiedPolynomials.substitute(minimalPolynomialOverUnramified,
						Collections.singletonList(unramifiedPolynomials.getPolynomial(root, unramified.one())));
			} else if (constantValuation.equals(Value.ONE) || minimalPolynomialOverUnramified.degree() == 1) {
				continue;
			} else {
				int ramification = minimalPolynomialOverUnramified.degree();
				int valuation = constantValuation.value();
				int gcd = BigInteger.valueOf(valuation).gcd(BigInteger.valueOf(ramification)).intValueExact();
				ramification /= gcd;
				valuation /= gcd;
				int alphaPower = BigInteger.valueOf(valuation).modInverse(BigInteger.valueOf(ramification))
						.intValueExact();
				int uniformizerPower = (alphaPower * valuation - 1) / ramification;
				alpha = divide(power(alpha, alphaPower), power(getEmbedding(base.uniformizer()), uniformizerPower));
				minimalPolynomialOverUnramified = maximalUnramifiedExtension.minimalPolynomial(alpha);
			}
		}
	}

	@Override
	public LocalRingExtension<B, Ext<B>, CompleteLocalFieldExtension<B, R, RE, RFE>, R, RE, RFE> ringOfIntegers() {
		return localRing;
	}

	@Override
	public RFE reduction() {
		return reduction;
	}

	@Override
	public RE reduce(Ext<B> t) {
		Vector<B> reductionGeneratorUniformizer = matrixAlgebra()
				.multiply(alphaToReductionGeneratorUniformizerBaseChange, asVector(t));
		List<R> reduced = new ArrayList<>();
		for (int i = 0; i < residueDegree; i++) {
			reduced.add(base.reduce(reductionGeneratorUniformizer.get(i + 1)));
		}
		return reduction.fromVector(new Vector<>(reduced));
	}
	
	@Override
	public Ext<B> upToUniformizerPower(Ext<B> t) {
		if (t.equals(zero())) {
			return zero();
		}
		Value value = valuation(t);
		Ext<B> uniformizerPower = power(uniformizer(), value.value());
		return divide(t, uniformizerPower);
	}

	@SuppressWarnings("unchecked")
	@Override
	public Ext<B> lift(RE s) {
		UnivariatePolynomialRing<Ext<B>> ring = getUnivariatePolynomialRing();
		return ring.evaluate(
				ring.getEmbedding(base.ringOfIntegers().liftUnivariatePolynomial(s.asPolynomial()), getEmbeddingMap()),
				reductionGenerator);
	}

	@Override
	public Value valuation(Ext<B> t) {
		if (alphaToReductionGeneratorUniformizerBaseChange == null) {
			return highAccuracyBase.valuation(highAccuracyAlgebra.determinant(asMatrix(t)));
		}
		Value value = Value.INFINITY;
		ValueGroup g = ValueGroup.g();
		Vector<B> reductionGeneratorUniformizer = matrixAlgebra()
				.multiply(alphaToReductionGeneratorUniformizerBaseChange, asVector(t));
		for (int i = 0; i < ramificationIndex; i++) {
			for (int j = 0; j < residueDegree; j++) {
				Value val = base.valuation(reductionGeneratorUniformizer.get(i * residueDegree + j + 1));
				val = g.multiply(ramificationIndex, val);
				val = g.operate(val, new Value(i));
				value = g.min(value, val);
			}
		}
		return value;
	}

	@Override
	public Real value(Ext<B> t) {
		Map<Real, Integer> roots = r.roots(base.value(norm(t)), degree());
		return r.abs(roots.keySet().iterator().next());
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	public int ramificationIndex() {
		return ramificationIndex;
	}

	public int residueDegree() {
		return residueDegree;
	}

	@Override
	public Ext<B> negative(Ext<B> t, int accuracy) {
		if (accuracy < this.getAccuracy()) {
			return round(negative(t), accuracy);
		} else if (accuracy == this.getAccuracy()) {
			return negative(t);
		}
		return withAccuracy(accuracy).negative(t);
	}

	@Override
	public Ext<B> inverse(Ext<B> t) {
		UnivariatePolynomialRing.ExtendedResultantResult<B> eres = highAccuracyPolynomials
				.extendedResultant(minimalPolynomial(), t.e);
		if (eres.getResultant().equals(base.zero())) {
			throw new ArithmeticException("Not invertible!");
		}
		return fromPolynomial(base.ringOfIntegers().roundUnivariatePolynomial(
				highAccuracyPolynomials.divideScalar(eres.getCoeff2(), eres.getGcd().leadingCoefficient()),
				base.getAccuracy()));
	}

	@Override
	public Ext<B> inverse(Ext<B> t, int accuracy) {
		if (accuracy < this.getAccuracy()) {
			return round(inverse(t), accuracy);
		} else if (accuracy == this.getAccuracy()) {
			return inverse(t);
		}
		return withAccuracy(accuracy).inverse(t);
	}

	@Override
	public Ext<B> divide(Ext<B> dividend, Ext<B> divisor) {
		if (dividend.equals(zero())) {
			return zero();
		}
		if (divisor.equals(zero())) {
			throw new ArithmeticException("Division by zero!");
		}
		return multiply(dividend, inverse(divisor,
				getAccuracy() + Math.max(0, -valuation(dividend).value() + valuation(divisor).value())));
	}

	@Override
	public PivotStrategy<Ext<B>> preferredPivotStrategy() {
		return new ValuationPivotStrategy<>(this.valuation());
	}

	@Override
	public boolean isInteger(Ext<B> t) {
		Vector<B> reductionGeneratorUniformizer = matrixAlgebra()
				.multiply(alphaToReductionGeneratorUniformizerBaseChange, asVector(t));
		for (int i = 0; i < ramificationIndex; i++) {
			for (int j = 0; j < residueDegree; j++) {
				if (!base.isInteger(reductionGeneratorUniformizer.get(i * residueDegree + j + 1))) {
					return false;
				}
			}
		}
		return true;
	}

	@Override
	public Ext<B> uniformizer() {
		return uniformizer;
	}

	@Override
	public List<Ext<B>> conjugates(Ext<B> s) {
		if (conjugates.containsKey(s)) {
			return conjugates.get(s);
		}
		UnivariatePolynomial<Ext<B>> mipo = getUnivariatePolynomialRing().getEmbedding(minimalPolynomial(s),
				getEmbeddingMap());
		if (ramificationIndex == 1) {
			List<Ext<B>> conjugates = new ArrayList<>();
			for (RE c : reduction().conjugates(reduce(s))) {
				conjugates.add(localRing.henselLift(mipo, c));
			}
			this.conjugates.put(s, conjugates);
			return conjugates;
		}
		List<Ext<B>> conjugates = new ArrayList<>();
		conjugates.addAll(roots(mipo).keySet());
		this.conjugates.put(s, conjugates);
		return conjugates(s);
	}

	@Override
	public GaloisGroup<B, Ext<B>, CompleteLocalFieldExtension<B, R, RE, RFE>> galoisGroup() {
		if (ramificationIndex == 1) {
			List<FieldAutomorphism<B, Ext<B>, CompleteLocalFieldExtension<B, R, RE, RFE>>> elements = new ArrayList<>();
			for (FieldAutomorphism<R, RE, RFE> reduced : reduction().galoisGroup()) {
				elements.add(new FieldAutomorphism<>(this, reduced.asArray()));
			}
			return new GaloisGroup<>(this, elements);
		}
		throw new UnsupportedOperationException("Not implemented");
	}

	@Override
	public Ext<B> round(Ext<B> t, int accuracy) {
		Vector<B> asVector = matrixAlgebra().multiply(alphaToReductionGeneratorUniformizerBaseChange, asVector(t));
		Ext<B> rounded = zero();
		Ext<B> uniformizerPower = one();
		for (int i = 0; i < ramificationIndex; i++) {
			Ext<B> reductionGeneratorPower = one();
			for (int j = 0; j < residueDegree; j++) {
				Ext<B> generator = multiply(uniformizerPower, reductionGeneratorPower);
				rounded = add(rounded, multiply(getEmbedding(base.round(asVector.get(residueDegree * i + j + 1),
						(accuracy - i + ramificationIndex - 1) / ramificationIndex)), generator));
				reductionGeneratorPower = multiply(reductionGeneratorPower, reductionGenerator);
			}
			uniformizerPower = multiply(uniformizerPower, uniformizer);
		}
		return rounded;
	}

	@Override
	public int getAccuracy() {
		return base.getAccuracy() * ramificationIndex;
	}

	@Override
	public CompleteLocalFieldExtension<B, R, RE, RFE> withAccuracy(int accuracy) {
		return new CompleteLocalFieldExtension<>(minimalPolynomial(), base.withAccuracy(accuracy), reduction);
	}
	

	@Override
	public OtherVersion<Ext<B>, ?, RE, ?> exact() {
		return null;
	}
		
	@Override
	public OtherVersion<Ext<B>, Ext<B>, RE, CompleteLocalFieldExtension<B, R, RE, RFE>> complete(int accuracy) {
		return new OtherVersion<>(withAccuracy(accuracy), new Identity<>(), new Identity<>());
	}


	@Override
	public Ext<B> getRandomInteger() {
		Ext<B> result = zero();
		Ext<B> uniformizerPower = one();
		for (int i = 0; i < ramificationIndex; i++) {
			Ext<B> reductionGeneratorPower = one();
			for (int j = 0; j < residueDegree; j++) {
				result = add(result,
						multiply(getEmbedding(base.getRandomInteger()), uniformizerPower, reductionGeneratorPower));
				reductionGeneratorPower = multiply(reductionGeneratorPower, reductionGenerator);
			}
			uniformizerPower = multiply(uniformizerPower, uniformizer);
		}
		return result;
	}
	
	@Override
	public FactorizationResult<Polynomial<Ext<B>>, Ext<B>> factorization(UnivariatePolynomial<Ext<B>> t) {
		return ringOfIntegers().factorization(t, true, getAccuracy());
	}

	@Override
	protected Ext<B> fromSmallDegreePolynomial(UnivariatePolynomial<B> polynomial) {
		return new Ext<>(polynomial);
	}

	@Override
	public CompleteLocalFieldExtension<B, R, RE, RFE> makeExtension(UnivariatePolynomial<B> minimalPolynomial) {
		return new CompleteLocalFieldExtension<>(minimalPolynomial, base, reduction);
	}

	@Override
	protected CompleteLocalFieldExtension<B, R, RE, RFE> asExtensionType() {
		return this;
	}
}
