package fields.local;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.interfaces.Algebra;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DedekindRingExtension;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationFieldExtension;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public class LocalRingExtension<B extends Element<B>, S extends Element<S>, E extends AlgebraicExtensionElement<B, E>, FE extends DiscreteValuationFieldExtension<B, S, E,FE, R, RE, RFE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>>
		extends LocalRingImplementation<E, RE> implements Algebra<B, E>, DiscreteValuationRing<E, RE>, DedekindRingExtension<B, B, S, E, E, R, RE, RFE, FE, FE> {
	private DiscreteValuationField<B, S> base;
	private FE extension;
	private RFE reductionExtension;
	private List<E> integralBasis;
	private Matrix<B> fromPowerToIntegralBaseChange;
	private FreeModule<B> asFreeModule;

	public LocalRingExtension(FE extension, DiscreteValuationField<B, S> base,  RFE reductionExtension,
			List<E> integralBasis, Matrix<B> fromPowerToIntegralBaseChange) {
		super(extension, base.ringOfIntegers().toString() + "[X]/(" + extension.minimalPolynomial() + ")");
		this.base = base;
		this.extension = extension;
		this.reductionExtension = reductionExtension;
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
	public List<Vector<B>> nonTrivialCombinations(List<E> s) {
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
	public List<Vector<B>> getSyzygies() {
		return Collections.emptyList();
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
	
	@Override
	public FE localField() {
		return extension;
	}
	
	@Override
	public FE quotientField() {
		return localField();
	}
	
	@Override
	public LocalRingExtension<B, S, E, FE, R, RE, RFE> localize(Ideal<E> maximalIdeal) {
		return this;
	}
	
	@Override
	public FE localizeAndQuotient(Ideal<E> maximalIdeal) {
		return localField();
	}

	private Matrix<B> asMatrix(List<E> m) {
		List<Vector<B>> result = new ArrayList<>();
		for (E s : m) {
			result.add(asVector(s));
		}
		return Matrix.fromColumns(result);
	}
	
	@Override
	public boolean isSubModuleMemberModule(List<E> m, E b) {
		Matrix<B> matrix = asMatrix(m);
		return isSubModuleMemberModule(matrix.getModule(getRing()), matrix, asVector(b));
	}
	
	@Override
	public boolean isSubModuleMemberModule(MatrixModule<B> module, Matrix<B> m, Vector<B> b) {
		return getRing().isSubModuleMember(module, m, b);
	}
	
	@Override
	public E asSubModuleMemberModule(List<E> m, E b) {
		Matrix<B> matrix = asMatrix(m);
		return fromVector(asSubModuleMemberModule(matrix.getModule(getRing()), matrix, asVector(b)));
	}
	
	@Override
	public Vector<B> asSubModuleMemberModule(MatrixModule<B> module, Matrix<B> m, Vector<B> b) {
		return getRing().asSubModuleMember(module, m, b);
	}
	
	@Override
	public List<Vector<B>> syzygyProblemModule(List<E> m) {
		Matrix<B> matrix = asMatrix(m);
		return syzygyProblemModule(matrix.getModule(getRing()), matrix);
	}
	
	@Override
	public List<Vector<B>> syzygyProblemModule(MatrixModule<B> module, Matrix<B> m) {
		return getRing().syzygyProblem(module, m);
	}
	
	@Override
	public List<E> simplifySubModuleGeneratorsModule(List<E> m) {
		Matrix<B> matrix = asMatrix(m);
		List<Vector<B>> simplified = simplifySubModuleGeneratorsModule(matrix.getModule(getRing()), matrix);
		List<E> result = new ArrayList<>();
		for (Vector<B> vector : simplified) {
			result.add(fromVector(vector));
		}
		return result;
	}
	
	@Override
	public List<Vector<B>> simplifySubModuleGeneratorsModule(MatrixModule<B> module, Matrix<B> m) {
		return getRing().simplifySubModuleGenerators(module, m);
	}
}
