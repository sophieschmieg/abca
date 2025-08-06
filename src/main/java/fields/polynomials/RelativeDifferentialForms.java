package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.AsSubModuleMemberPreperation;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.DifferentialForms.DifferentialForm;
import fields.vectors.ExteriorProduct.WedgeProduct;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public class RelativeDifferentialForms<T extends Element<T>>
		extends AbstractModule<CoordinateRingElement<T>, RelativeDifferentialForms<T>.RelativeDifferentialForm> {
	public class RelativeDifferentialForm extends AbstractElement<RelativeDifferentialForm> {
		private DifferentialForm<T> differentialForm;

		private RelativeDifferentialForm(DifferentialForm<T> differentialForm) {
			this.differentialForm = differentialForm;
		}

		@Override
		public String toString() {
			return differentialForm.toString();
		}

		@Override
		public int compareTo(RelativeDifferentialForm o) {
			return differentialForm.compareTo(o.differentialForm);
		}
	}

	private CoordinateRing<T> base;
	private List<Polynomial<T>> map;
	private PolynomialRing<T> polynomialRing;
	private CoordinateRing<T> coordinateRing;
	private DifferentialForms<T> differentialForms;
	private Map<Integer, AsSubModuleMemberPreperation<CoordinateRingElement<T>, ?>> generatorPreparations;
	private FreeModule<CoordinateRingElement<T>> asFreeModule;

	public RelativeDifferentialForms(CoordinateRing<T> coordinateRing) {
		this(coordinateRing,
				coordinateRing.getPolynomialRing().getRing().getUnivariatePolynomialRing().getZeroIdeal().divideOut(),
				Collections.singletonList(coordinateRing.getPolynomialRing().zero()));
	}

	public RelativeDifferentialForms(CoordinateRing<T> coordinateRing, CoordinateRing<T> base,
			List<Polynomial<T>> map) {
		if (base.getPolynomialRing().numberOfVariables() != map.size()) {
			throw new ArithmeticException("Not enough maps for the variables");
		}
		this.base = base;
		this.map = map;
		this.coordinateRing = coordinateRing;
		this.polynomialRing = coordinateRing.getPolynomialRing();
		this.differentialForms = polynomialRing.differentialForms();
		this.asFreeModule = new FreeModule<>(coordinateRing, polynomialRing.numberOfVariables());
		this.generatorPreparations = new TreeMap<>();

	}

	private AsSubModuleMemberPreperation<CoordinateRingElement<T>, ?> generatorPreparation(int degree) {
		if (!generatorPreparations.containsKey(degree)) {
			List<Polynomial<T>> generatorPolynomials = new ArrayList<>();
			for (Polynomial<T> generator : coordinateRing.getIdeal().generators()) {
				generatorPolynomials.add(generator);
			}
			for (int i = 0; i < base.getPolynomialRing().numberOfVariables(); i++) {
				CoordinateRingElement<T> additionalGenerator = coordinateRing.zero();
				CoordinateRingElement<T> var = base.getVar(i + 1);
				for (Monomial m : var.getElement().monomials()) {
					CoordinateRingElement<T> term = coordinateRing.getEmbedding(var.getElement().coefficient(m));
					for (int j = 0; j < base.getPolynomialRing().numberOfVariables(); j++) {
						term = coordinateRing.multiply(
								coordinateRing.power(coordinateRing.getEmbedding(map.get(j)), m.exponents()[j]), term);
					}
					additionalGenerator = coordinateRing.add(term, additionalGenerator);
				}
				generatorPolynomials.add(additionalGenerator.getElement());
			}
			List<Vector<CoordinateRingElement<T>>> generatorVectors = new ArrayList<>();
			for (Polynomial<T> generator : generatorPolynomials) {
				for (DifferentialForm<T> diff : differentialForms.getGradedModuleGenerators(degree - 1)) {
					DifferentialForm<T> differentialForm = differentialForms
							.derivative(differentialForms.scalarMultiply(generator, diff));
					generatorVectors.add(Vector.mapVector(coordinateRing.getEmbeddingMap(),
							differentialForms.asGradedVector(differentialForm, degree)));
				}
			}

			Matrix<CoordinateRingElement<T>> generators = Matrix.fromColumns(generatorVectors);
			MatrixModule<CoordinateRingElement<T>> generatorModule = generators.getModule(coordinateRing);
			this.generatorPreparations.put(degree,
					coordinateRing.prepareAsSubModuleMember(generatorModule, generators));

		}
		return generatorPreparations.get(degree);
	}

	public RelativeDifferentialForm getElement(DifferentialForm<T> differentialForm) {
		Map<WedgeProduct, Polynomial<T>> map = differentialForms.asMap(differentialForm);
		Map<WedgeProduct, Polynomial<T>> reducedMap = new TreeMap<>();
		int power = -1;
		for (WedgeProduct product : map.keySet()) {
			if (power == -1) {
				power = product.getVariables().size();
			} else if (power != product.getVariables().size()) {
				throw new ArithmeticException("Mixed differential!");
			}
			reducedMap.put(product, coordinateRing.getEmbedding(map.get(product)).getElement());
		}
		if (power == -1) {
			return new RelativeDifferentialForm(differentialForms.zero());
		}
		Vector<Polynomial<T>> vector = differentialForms.asGradedVector(differentialForms.fromMap(reducedMap), power);
		DifferentialForm<T> reduced = differentialForms.fromGradedVector(
				Vector.mapVector(coordinateRing.getAsElementMap(), coordinateRing.residueModuloSubmodule(
						generatorPreparation(power), Vector.mapVector(coordinateRing.getEmbeddingMap(), vector))),
				power);
		return new RelativeDifferentialForm(reduced);
	}

	public DifferentialForm<T> asDifferentialForm(RelativeDifferentialForm t) {
		return t.differentialForm;
	}

	public RelativeDifferentialForm derivative(CoordinateRingElement<T> t) {
		return getElement(differentialForms.derivative(t.getElement()));
	}

	public RelativeDifferentialForm derivative(RelativeDifferentialForm t) {
		return getElement(differentialForms.derivative(asDifferentialForm(t)));
	}

	@Override
	public CoordinateRing<T> getRing() {
		return coordinateRing;
	}

	@Override
	public RelativeDifferentialForm zero() {
		return fromVector(asFreeModule.zero());
	}

	@Override
	public RelativeDifferentialForm add(RelativeDifferentialForm s1, RelativeDifferentialForm s2) {
		return fromVector(asFreeModule.add(asVector(s1), asVector(s2)));
	}

	@Override
	public RelativeDifferentialForm negative(RelativeDifferentialForm s) {
		return fromVector(asFreeModule.negative(asVector(s)));
	}

	@Override
	public RelativeDifferentialForm scalarMultiply(CoordinateRingElement<T> t, RelativeDifferentialForm s) {
		return fromVector(asFreeModule.scalarMultiply(t, asVector(s)));
	}

	@Override
	public boolean isFree() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Ideal<CoordinateRingElement<T>> annihilator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isLinearIndependent(List<RelativeDifferentialForm> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isGeneratingModule(List<RelativeDifferentialForm> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<Vector<CoordinateRingElement<T>>> nonTrivialCombinations(List<RelativeDifferentialForm> s) {
		// TODO Auto-generated method stub
		return null;
	}

	public List<RelativeDifferentialForm> getGradedModuleGenerators(int degree) {
		List<RelativeDifferentialForm> result = new ArrayList<>();
		for (DifferentialForm<T> unitVector : differentialForms.getGradedModuleGenerators(degree)) {
			RelativeDifferentialForm generator = getElement(unitVector);
			if (!generator.equals(zero())) {
				result.add(generator);
			}
		}
		return result;
	}

	public List<Vector<CoordinateRingElement<T>>> getSyzygies(int degree) {
		return generatorPreparation(degree).getMatrix().asColumnList();
	}

	@Override
	public List<RelativeDifferentialForm> getModuleGenerators() {
		List<RelativeDifferentialForm> result = new ArrayList<>();
		for (DifferentialForm<T> unitVector : differentialForms.getModuleGenerators()) {
			RelativeDifferentialForm generator = getElement(unitVector);
			if (!generator.equals(zero())) {
				result.add(generator);
			}
		}
		return result;
	}
	
	@Override
	public List<Vector<CoordinateRingElement<T>>> getSyzygies() {
		// TODO Auto-generated method stub
		return null;
	}

	public RelativeDifferentialForm fromPolynomialVector(Vector<Polynomial<T>> vector, int degree) {
		return getElement(differentialForms.fromGradedVector(vector, degree));
	}
	
	public Vector<CoordinateRingElement<T>> asGradedVector(RelativeDifferentialForm s, int degree) {
		return Vector.mapVector(coordinateRing.getEmbeddingMap(),  differentialForms.asGradedVector(asDifferentialForm(s), degree));
	}
	
	@Override
	public Vector<CoordinateRingElement<T>> asVector(RelativeDifferentialForms<T>.RelativeDifferentialForm s) {
		return Vector.mapVector(coordinateRing.getEmbeddingMap(),  differentialForms.asVector(asDifferentialForm(s)));
	}

	@Override
	public Exactness exactness() {
		return coordinateRing.exactness();
	}

	@Override
	public RelativeDifferentialForm getRandomElement() {
		return fromVector(asFreeModule.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return coordinateRing.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return BigInteger.valueOf(-1);
	}

	@Override
	public Iterator<RelativeDifferentialForm> iterator() {
		throw new UnsupportedOperationException();
	}

}
