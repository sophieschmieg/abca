package fields.local;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import fields.interfaces.Algebra;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;

public class LocalRingExtension<B extends Element<B>, S extends Element<S>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE> & DiscreteValuationField<E, RE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>>
		extends LocalRingImplementation<E, RE> implements Algebra<B, E>, DiscreteValuationRing<E, RE> {
	private DiscreteValuationField<B, S> base;
	private FE extension;
	private RFE reductionExtension;
	private List<E> integralBasis;
	private Matrix<B> fromPowerToIntegralBaseChange;
	private FreeModule<B> asFreeModule;

	public LocalRingExtension(FE extension, DiscreteValuationField<B, S> base, OkutsuType<B, S, R, RE, RFE> type,
			List<E> integralBasis, Matrix<B> fromPowerToIntegralBaseChange) {
		super(extension, base.ringOfIntegers().toString() + "[X]/(" + extension.minimalPolynomial() + ")");
		this.base = base;
		this.extension = extension;
		this.reductionExtension = type.reduction().extension();
		this.integralBasis = integralBasis;
		this.fromPowerToIntegralBaseChange = fromPowerToIntegralBaseChange;
		this.asFreeModule = new FreeModule<>(base, extension.degree());
		if (this.extension.residueField() != this.reductionExtension) {
			throw new ArithmeticException("Not the reduction of the extension!");
		}
		if (this.extension.getBaseField() != this.base) {
			throw new ArithmeticException("Not the base field of the extension!");
		}

	}

	@Override
	public DiscreteValuationRing<B, S> getRing() {
		return base.ringOfIntegers();
	}

	@Override
	public RFE reduction() {
		return reductionExtension;
	}

	@Override
	public E scalarMultiply(B t, E s) {
		return extension.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}
	
	@Override
	public boolean isTorsionFree() {
		return true;
	}

	@Override
	public Ideal<B> annihilator() {
		return base.getZeroIdeal();
	}

	@Override
	public E scalarMultiply(int n, E s) {
		return extension.scalarMultiply(n, s);
	}

	@Override
	public E scalarMultiply(BigInteger n, E s) {
		return extension.scalarMultiply(n, s);
	}

	@Override
	public E scalarMultiply(B t1, B t2, E s) {
		return extension.scalarMultiply(t1, t2, s);
	}

	@Override
	public E scalarMultiply(int n, B t, E s) {
		return extension.scalarMultiply(n, t, s);
	}

	@Override
	public E scalarMultiply(BigInteger n, B t, E s) {
		return extension.scalarMultiply(n, t, s);
	}

	@Override
	public E scalarMultiply(int n, B t1, B t2, E s) {
		return extension.scalarMultiply(n, t1, t2, s);
	}

	@Override
	public E scalarMultiply(BigInteger n, B t1, B t2, E s) {
		return extension.scalarMultiply(n, t1, t2, s);
	}

	@Override
	public boolean isLinearIndependent(List<E> s) {
		List<Vector<B>> asVectors = new ArrayList<>();
		for (E e : s) {
			asVectors.add(asVector(e));
		}
		return asFreeModule.isLinearIndependent(asVectors);
	}

	@Override
	public boolean isGeneratingModule(List<E> s) {
		List<Vector<B>> asVectors = new ArrayList<>();
		for (E e : s) {
			asVectors.add(asVector(e));
		}
		return asFreeModule.isGeneratingModule(asVectors);
	}

	@Override
	public List<B> nonTrivialCombination(List<E> s) {
		return nonTrivialCombinations(s).get(0);
	}

	@Override
	public List<List<B>> nonTrivialCombinations(List<E> s) {
		List<Vector<B>> asVectors = new ArrayList<>();
		for (E e : s) {
			asVectors.add(asVector(e));
		}
		return asFreeModule.nonTrivialCombinations(asVectors);
	}

	@Override
	public List<E> getModuleGenerators() {
		return integralBasis;
	}

	@Override
	public Vector<B> asVector(E s) {
		return extension.matrixAlgebra().multiply(fromPowerToIntegralBaseChange, extension.asVector(s));
	}

	@Override
	public E fromVector(Vector<B> vector) {
		E result = zero();
		List<E> generators = getModuleGenerators();
		for (int i = 0; i < generators.size(); i++) {
			result = add(result, scalarMultiply(vector.get(i + 1), generators.get(i)));
		}
		return result;
	}

	@Override
	public E getEmbedding(B t) {
		return extension.getEmbedding(t);
	}

	@Override
	public boolean isGeneratingAlgebra(List<E> s) {
		List<E> includingPowers = new ArrayList<>();
		includingPowers.add(one());
		for (E e : s) {
			E power = one();
			for (int i = 1; i < extension.degree(); i++) {
				power = multiply(power, e);
				includingPowers.add(power);
			}
		}
		return isGeneratingModule(includingPowers);
	}

	@Override
	public List<E> getAlgebraGenerators() {
		return getModuleGenerators();
	}
}
