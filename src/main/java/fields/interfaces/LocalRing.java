package fields.interfaces;

import java.util.List;

import fields.helper.FieldEmbedding;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Field.Extension;
import fields.interfaces.LocalField.Valuation;
import fields.local.Value;
import util.Identity;
import util.Pair;

public interface LocalRing<T extends Element<T>, S extends Element<S>> extends ValueRing<T>, DedekindRing<T, T, S> {
	boolean isComplete();

	public boolean isElement(T t);

	public Value valuation(T t);

	public Value valuationOfUnivariatePolynomial(Polynomial<T> t);

	@Override
	public default Value valuation(T t, Ideal<T> maximalIdeal) {
		return valuation(t);
	}

	public Valuation<T> valuation();

	public T uniformizer();

	public LocalField<T, S> fieldOfFractions();

	public default MathMap<T, T> embedding() {
		return new Identity<>();
	}

	public default boolean isInteger(T t) {
		return isElement(t);
	}

	public default T asInteger(T t) {
		if (!isElement(t)) {
			throw new ArithmeticException("Not an element");
		}
		return t;
	}

	public Field<S> reduction();

	@Override
	public default Field<S> reduction(Ideal<T> maximalIdeal) {
		return reduction();
	}

	@Override
	public default LocalRing<T, S> localize(Ideal<T> maximalIdeal) {
		return this;
	}

	public S reduce(T t);

	@Override
	public default S reduce(T t, Ideal<T> maximalIdeal) {
		return reduce(t);
	}

	public Polynomial<S> reducePolynomial(Polynomial<T> t);

	public T upToUniformizerPower(T t);

	public UnivariatePolynomial<S> reduceUnivariatePolynomial(UnivariatePolynomial<T> t);

	public T lift(S s);

	@Override
	public default T lift(S s, Ideal<T> maximalIdeal) {
		return lift(s);
	}

	public Polynomial<T> liftPolynomial(Polynomial<S> t);

	public UnivariatePolynomial<T> liftUnivariatePolynomial(Polynomial<S> t);

	public T henselLift(UnivariatePolynomial<T> f, S aReduced);

	public interface OkutsuType<T extends Element<T>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> {
		public UnivariatePolynomial<T> getPolynomial();

		public UnivariatePolynomial<T> representative();

		public UnivariatePolynomial<T> phi();

		public Fraction lambda();

		public UnivariatePolynomial<RE> psi();

		public int level();

		public OkutsuType<T, S, R, RE, RFE> previousLevel();

		public boolean complete();

		public boolean leaf();

		public OkutsuType<T, S, R, RE, RFE> updateLeaf(UnivariatePolynomial<T> representative);

		public Value value();

		public Value quality();

		public Value precision();

		public int exponent();

		public int ramificationIndex();

		public int residueDegree();

		public int[] uniformizer();

		public Extension<S, R, RE, RFE> reduction();

		public FieldEmbedding<R, RE, RFE> embeddedReductionExtension();

		public Value valuation(UnivariatePolynomial<T> t);

		public RE reduce(UnivariatePolynomial<T> t);

		public UnivariatePolynomial<RE> reduceAsPolynomial(UnivariatePolynomial<T> t);

		public UnivariatePolynomial<T> lift(RE t, int value);

		public UnivariatePolynomial<T> liftAsPolynomial(UnivariatePolynomial<RE> t, int value);

		public boolean divides(UnivariatePolynomial<T> t);

		public UnivariatePolynomial<T> divisorPolynomial(int m);
	}

	public class TheMontesResult<T extends Element<T>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> {
		private List<OkutsuType<T, S, R, RE, RFE>> types;
		private UnivariatePolynomial<T> f;

		public TheMontesResult(List<OkutsuType<T, S, R, RE, RFE>> types, UnivariatePolynomial<T> f) {
			this.types = types;
			this.f = f;
		}

		public List<OkutsuType<T, S, R, RE, RFE>> getTypes() {
			return types;
		}

		public UnivariatePolynomial<T> getPolynomial() {
			return f;
		}
		
		@Override
		public String toString() {
			return types.toString();
		}
	}

	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> TheMontesResult<T, S, R, RE, RFE> theMontesAlgorithm(
			UnivariatePolynomial<T> t, Extension<S, R, RE, RFE> reductionAsExtension);

	public TheMontesResult<T, S, ?, ?, ?> theMontesAlgorithm(UnivariatePolynomial<T> f);

	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> UnivariatePolynomial<T> invertInType(
			UnivariatePolynomial<T> t, OkutsuType<T, S, R, RE, RFE> okutsuType, int precision);

	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> OkutsuType<T, S, R, RE, RFE> singleFactorLifting(
			OkutsuType<T, S, R, RE, RFE> type, int accuracy);

	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> List<List<UnivariatePolynomial<T>>> integralBasis(
			Polynomial<T> minimalPolynomial, TheMontesResult<T, S, R, RE, RFE> theMontes, boolean reduced);

	public List<UnivariatePolynomial<T>> triagonalizeIntegralBasis(Polynomial<T> minimalPolynomial,
			List<List<UnivariatePolynomial<T>>> integralBasis);

	public static class HenselLiftFactorResult<T extends Element<T>> {
		private UnivariatePolynomial<T> liftOfFactor;
		private UnivariatePolynomial<T> liftOfCofactor;
		private UnivariatePolynomial<T> factorCoefficient;
		private UnivariatePolynomial<T> cofactorCoefficient;

		public HenselLiftFactorResult(UnivariatePolynomial<T> liftOfFactor, UnivariatePolynomial<T> liftOfCofactor,
				UnivariatePolynomial<T> factorCoefficient, UnivariatePolynomial<T> cofactorCoefficient) {
			this.liftOfFactor = liftOfFactor;
			this.liftOfCofactor = liftOfCofactor;
			this.factorCoefficient = factorCoefficient;
			this.cofactorCoefficient = cofactorCoefficient;
		}

		public UnivariatePolynomial<T> getLiftOfFactor() {
			return liftOfFactor;
		}

		public UnivariatePolynomial<T> getLiftOfCofactor() {
			return liftOfCofactor;
		}

		public UnivariatePolynomial<T> getFactorCoefficient() {
			return factorCoefficient;
		}

		public UnivariatePolynomial<T> getCofactorCoefficient() {
			return cofactorCoefficient;
		}
	}

	public UnivariatePolynomial<T> henselLiftFactor(UnivariatePolynomial<T> f, UnivariatePolynomial<S> gReduced);

	public UnivariatePolynomial<T> henselLiftFactor(UnivariatePolynomial<T> f, UnivariatePolynomial<S> gReduced,
			int accuracy);

	public HenselLiftFactorResult<T> extendedHenselLiftFactor(UnivariatePolynomial<T> f,
			UnivariatePolynomial<S> gReduced);

	public HenselLiftFactorResult<T> extendedHenselLiftFactor(UnivariatePolynomial<T> f,
			UnivariatePolynomial<S> gReduced, int accuracy);

	public T henselLift(UnivariatePolynomial<T> f, S aReduced, int accuracy);

	public T henselLiftWithInitialLift(UnivariatePolynomial<T> f, T initialLift);

	public T henselLiftWithInitialLift(UnivariatePolynomial<T> f, T initialLift, int accuracy);

	public FactorizationResult<Polynomial<T>, T> henselLiftFactorization(UnivariatePolynomial<T> f, int accuracy);

	public T round(T t, int accuracy);

	public Polynomial<T> roundPolynomial(Polynomial<T> t, int accuracy);

	public UnivariatePolynomial<T> roundUnivariatePolynomial(UnivariatePolynomial<T> t, int accuracy);

	public boolean hasGoodReduction(UnivariatePolynomial<T> t);

	public boolean hasIrreducibleGoodReduction(UnivariatePolynomial<T> t);

	public boolean isIntegral(UnivariatePolynomial<T> t);

	public boolean isEisenstein(UnivariatePolynomial<T> t);

	public UnivariatePolynomial<T> integralMinimalPolynomial(UnivariatePolynomial<T> minimalPolynomial);

	public Pair<UnivariatePolynomial<T>, UnivariatePolynomial<T>> integralPolynomial(UnivariatePolynomial<T> t);

	public List<UnivariatePolynomial<T>> localExtensions(UnivariatePolynomial<T> minimalPolynomial);

	public boolean isIrreducible(UnivariatePolynomial<T> t);

	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t, boolean forField, int accuracy);

	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t, int accuracy);

	Module<T, T> fieldOfFractionsAsModule();

}
