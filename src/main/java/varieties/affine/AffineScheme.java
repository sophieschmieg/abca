package varieties.affine;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.IdealResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.DifferentialForms;
import fields.polynomials.FieldExtensionUnivariatePolynomialRing;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.FunctionMathMap;
import varieties.AbstractMorphism;
import varieties.AbstractScheme;
import varieties.GeneralRationalFunction;
import varieties.Morphism;
import varieties.SimpleFunctionField;
import varieties.SimpleFunctionField.SimpleFunctionFieldFromCoordinateRing;
import varieties.SimpleFunctionField.SimpleRationalFunction;
import varieties.SpectrumOfField;
import varieties.SpectrumOfField.SingletonPoint;
import varieties.affine.AffineMorphism.OpenImmersionResult;

public class AffineScheme<T extends Element<T>> extends AbstractScheme<T, AffinePoint<T>> {
	private CoordinateRing<T> coordinateRing;
	private Field<T> field;
	private AffineCover<T> cover;
	private Map<AffinePoint<T>, LocalizedCoordinateRing<T>> localRings;
	private AffineSimplification<T> simplification;

	public AffineScheme(Field<T> field, CoordinateRing<T> coordinateRing) {
		if (!coordinateRing.getPolynomialRing().getRing().equals(field)) {
			throw new ArithmeticException("Incompatible parameters");
		}
		this.coordinateRing = coordinateRing;
		this.field = field;
		this.cover = new AffineCover<>(Collections.singletonList(this),
				Collections.singletonList(Collections.singletonList(identityMorphism())), false);
		this.localRings = new TreeMap<>();
	}

	public static class IntersectionResult<T extends Element<T>> {
		private AffineScheme<T> intersection;
		private AffineMorphism<T> firstEmbedding;
		private AffineMorphism<T> secondEmbedding;

		public AffineScheme<T> getIntersection() {
			return intersection;
		}

		public AffineMorphism<T> getFirstEmbedding() {
			return firstEmbedding;
		}

		public AffineMorphism<T> getSecondEmbedding() {
			return secondEmbedding;
		}
	}

	public static <T extends Element<T>> IntersectionResult<T> intersect(AffineScheme<T> t1, AffineScheme<T> t2) {
		PolynomialRing<T> polynomialRing = t1.coordinateRing.getPolynomialRing();
		if (!polynomialRing.equals(t2.coordinateRing.getPolynomialRing())) {
			throw new ArithmeticException("different polynomial rings!");
		}
		PolynomialIdeal<T> ideal = polynomialRing.add(t1.coordinateRing.getIdeal(), t2.coordinateRing.getIdeal());
		List<Polynomial<T>> identity = new ArrayList<>();
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			identity.add(polynomialRing.getVar(i + 1));
		}
		IntersectionResult<T> result = new IntersectionResult<>();
		result.intersection = new AffineScheme<>(t1.getField(), ideal.divideOut());
		result.firstEmbedding = AffineMorphism.fromPolynomials(result.intersection, t1, identity);
		result.secondEmbedding = AffineMorphism.fromPolynomials(result.intersection, t2, identity);
		return result;
	}

//	public static <T extends Element<T>> AffineMorphism<T> restrictAwayFrom(AffineScheme<T> variety,
//			CoordinateIdeal<T> ideal) {
//		List<CoordinateRingElement<T>> constraints = ideal.generators();
//		int numVars = variety.coordinateRing.getPolynomialRing().numberOfVariables();
//		PolynomialRing<T> polynomialRing = AbstractPolynomialRing.getPolynomialRing(variety.field,
//				numVars + constraints.size(), variety.coordinateRing.getPolynomialRing().getComparator());
//		List<Polynomial<T>> idealGenerators = new ArrayList<>();
//		for (Polynomial<T> generator : variety.coordinateRing.getIdeal().generators()) {
//			idealGenerators.add(polynomialRing.getEmbedding(generator));
//		}
//		for (int i = 0; i < constraints.size(); i++) {
//			Polynomial<T> t = polynomialRing.getEmbedding(constraints.get(i).getElement());
//			idealGenerators.add(polynomialRing.subtract(
//					polynomialRing.multiply(t, polynomialRing.getVar(i + 1 + numVars)), polynomialRing.one()));
//		}
//		PolynomialIdeal<T> polynomialIdeal = polynomialRing.getIdeal(idealGenerators);
//		List<Polynomial<T>> identity = new ArrayList<>();
//		for (int i = 0; i < numVars; i++) {
//			identity.add(polynomialRing.getVar(i + 1));
//		}
//		AffineScheme<T> restricted = new AffineScheme<>(variety.getField(), polynomialIdeal.divideOut());
//		return AffineMorphism.fromPolynomials(restricted, variety, identity);
//	}

	public static class UnionResult<T extends Element<T>> {
		private AffineScheme<T> union;
		private AffineMorphism<T> firstEmbedding;
		private AffineMorphism<T> secondEmbedding;

		public AffineScheme<T> getUnion() {
			return union;
		}

		public AffineMorphism<T> getFirstEmbedding() {
			return firstEmbedding;
		}

		public AffineMorphism<T> getSecondEmbedding() {
			return secondEmbedding;
		}
	}

	public static <T extends Element<T>> UnionResult<T> union(AffineScheme<T> t1, AffineScheme<T> t2) {
		PolynomialRing<T> polynomialRing = t1.coordinateRing.getPolynomialRing();
		if (!polynomialRing.equals(t2.coordinateRing.getPolynomialRing())) {
			throw new ArithmeticException("different polynomial rings!");
		}
		PolynomialIdeal<T> ideal = polynomialRing.intersect(t1.coordinateRing.getIdeal(), t2.coordinateRing.getIdeal());
		List<Polynomial<T>> identity = new ArrayList<>();
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			identity.add(polynomialRing.getVar(i + 1));
		}
		UnionResult<T> result = new UnionResult<>();
		result.union = new AffineScheme<>(t1.getField(), ideal.divideOut());
		result.firstEmbedding = AffineMorphism.fromPolynomials(t1, result.union, identity);
		result.secondEmbedding = AffineMorphism.fromPolynomials(t2, result.union, identity);
		return result;
	}

	@Override
	public String toString() {
		return "Spec " + coordinateRing;
	}

	public CoordinateRing<T> getCoordinateRing() {
		return coordinateRing;
	}

	public LocalizedCoordinateRing<T> localizedCoordinateRing(AffinePoint<T> point) {
		if (!localRings.containsKey(point)) {
			CoordinateIdeal<T> ideal = coordinateRing.getIdeal(point.asIdeal(coordinateRing.getPolynomialRing()));
			localRings.put(point, new LocalizedCoordinateRing<>(field, coordinateRing, ideal));
		}
		return localRings.get(point);
	}

	public boolean hasRationalPoint(AffinePoint<T> point) {
		return point.asIdeal(coordinateRing.getPolynomialRing()).contains(coordinateRing.getIdeal());
	}

	public Field<T> getField() {
		return field;
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof AffineScheme<?>)) {
			return false;
		}
		@SuppressWarnings("unchecked")
		AffineScheme<T> other = (AffineScheme<T>) obj;
		return coordinateRing.equals(other.coordinateRing);
	}

	@Override
	public AffinePoint<T> getRandomElement() {
		List<Polynomial<T>> toSolve = new ArrayList<>();
		PolynomialRing<T> polynomialRing = coordinateRing.getPolynomialRing();
		for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
			if (coordinateRing.boundVariables().contains(i + 1)) {
				continue;
			}
			toSolve.add(polynomialRing.subtract(polynomialRing.getVar(i + 1),
					polynomialRing.getEmbedding(field.getRandomElement())));
		}
		toSolve.addAll(coordinateRing.getIdeal().generators());
		List<AffinePoint<T>> points = polynomialRing.solve(toSolve);
		return points.get(new Random().nextInt(points.size()));
	}

	@Override
	public boolean isFinite() {
		return field.isFinite() || dimension() == 0;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		if (dimension() == 0) {
			return BigInteger.valueOf(coordinateRing.getPolynomialRing().solve(coordinateRing.getIdeal()).size());
		}
		throw new UnsupportedOperationException();
	}

	@Override
	public Iterator<AffinePoint<T>> iterator() {
		if (dimension() < 0) {
			return Collections.emptyIterator();
		}
		PolynomialRing<T> polynomialRing = coordinateRing.getPolynomialRing();
		PolynomialRing<T> boundRing = AbstractPolynomialRing.getPolynomialRing(field,
				polynomialRing.numberOfVariables() - dimension(), polynomialRing.getComparator());
		if (dimension() == 0) {
			return polynomialRing.solve(coordinateRing.getIdeal()).iterator();
		}
		FreeModule<T> module = new FreeModule<>(field, dimension());
		int[] map = new int[polynomialRing.numberOfVariables()];
		int variable = 0;
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			if (coordinateRing.boundVariables().contains(i + 1)) {
				map[i] = variable;
				variable++;
				continue;
			}
			map[i] = -1;
		}
		return new Iterator<>() {
			Deque<AffinePoint<T>> queue = new LinkedList<>();
			Iterator<Vector<T>> it = module.iterator();

			private boolean fillQueue() {
				while (queue.size() <= 1) {
					if (!it.hasNext()) {
						return !queue.isEmpty();
					}
					Vector<T> next = it.next();
					List<T> eval = new ArrayList<>();
					int variable = 1;
					for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
						if (coordinateRing.boundVariables().contains(i + 1)) {
							eval.add(null);
							continue;
						}
						eval.add(next.get(variable));
						variable++;
					}
					List<Polynomial<T>> evaluated = new ArrayList<>();
					for (Polynomial<T> generator : coordinateRing.getIdeal().generators()) {
						evaluated.add(boundRing.getEmbedding(polynomialRing.partiallyEvaluate(generator, eval), map));
					}
					List<AffinePoint<T>> solutions = boundRing.solve(evaluated);
					for (AffinePoint<T> solution : solutions) {
						List<T> point = new ArrayList<>();
						int boundVariable = 1;
						int freeVariable = 1;
						for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
							if (coordinateRing.boundVariables().contains(i + 1)) {
								point.add(solution.getCoord(boundVariable));
								boundVariable++;
								continue;
							}
							point.add(next.get(freeVariable));
							freeVariable++;
						}
						queue.add(new AffinePoint<>(field, point));
					}
				}
				return true;
			}

			@Override
			public boolean hasNext() {
				return fillQueue();
			}

			@Override
			public AffinePoint<T> next() {
				if (!fillQueue()) {
					throw new NoSuchElementException("Empty iterator!");
				}
				return queue.poll();
			}
		};
	}

	@Override
	public AffineCover<T> getAffineCover() {
		return cover;
	}

	@Override
	public List<Integer> affineCoverIndex(AffinePoint<T> p) {
		return Collections.singletonList(0);
	}

	@Override
	public int recommendAffineCoverIndex(AffinePoint<T> p) {
		return 0;
	}

	@Override
	public AffineMorphism<T> identityMorphism() {
		List<Polynomial<T>> vars = new ArrayList<>();
		PolynomialRing<T> r = coordinateRing.getPolynomialRing();
		for (int i = 0; i < r.numberOfVariables(); i++) {
			vars.add(r.getVar(i + 1));
		}
		return AffineMorphism.fromPolynomials(this, this, vars);
	}

	@Override
	public Morphism<T, AffinePoint<T>, AffinePoint<T>> embedding(int coverIndex) {
		if (coverIndex != 0) {
			throw new IndexOutOfBoundsException();
		}
		return identityMorphism();
	}

	@Override
	public AffinePoint<T> asAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		if (affineCoverIndex != 0) {
			throw new IndexOutOfBoundsException();
		}
		return p;
	}

	@Override
	public AffinePoint<T> fromAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		if (affineCoverIndex != 0) {
			throw new IndexOutOfBoundsException();
		}
		return p;
	}

	@Override
	public Morphism<T, SingletonPoint, AffinePoint<T>> pointAsMorphism(AffinePoint<T> p) {
		SpectrumOfField<T> domain = new SpectrumOfField<>(field);
		return new AbstractMorphism<>() {

			@Override
			public AffinePoint<T> evaluate(SingletonPoint t) {
				return p;
			}

			@Override
			public SpectrumOfField<T> getDomain() {
				return domain;
			}

			@Override
			public AffineScheme<T> getRange() {
				return AffineScheme.this;
			}

			@Override
			public GeneralRationalFunction<T, SingletonPoint, AffinePoint<T>> restrict(SingletonPoint point) {
				AffineScheme<T> singleton = domain.getAffineCover().getCover().get(0);
				PolynomialRing<T> polynomials = singleton.coordinateRing.getPolynomialRing();
				List<Polynomial<T>> asPolynomials = new ArrayList<>();
				for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
					asPolynomials.add(polynomials.getEmbedding(p.getCoord(i + 1)));
				}
				return new GeneralRationalFunction<>(domain, AffineScheme.this,
						AffineMorphism.fromPolynomials(singleton, AffineScheme.this, asPolynomials), 0,
						singleton.identityMorphism(), 0);
			}

		};
	}

	@Override
	public List<AffineMorphism<T>> irreducibleComponents() {
		List<AffineMorphism<T>> result = new ArrayList<>();
		List<PolynomialIdeal<T>> components = coordinateRing.getPolynomialRing()
				.primaryDecomposition(coordinateRing.getIdeal()).getPrimaryIdeals();
		for (PolynomialIdeal<T> component : components) {
			AffineScheme<T> domain = new AffineScheme<>(field, component.divideOut());
			List<Polynomial<T>> polynomials = new ArrayList<>();
			for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
				polynomials.add(coordinateRing.getPolynomialRing().getVar(i + 1));
			}
			result.add(AffineMorphism.fromPolynomials(domain, this, polynomials));
		}
		return result;
	}

	@Override
	public AffineMorphism<T> reduced() {
		AffineScheme<T> domain = new AffineScheme<>(field,
				coordinateRing.getPolynomialRing().radical(coordinateRing.getIdeal()).divideOut());
		List<Polynomial<T>> polynomials = new ArrayList<>();
		for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
			polynomials.add(coordinateRing.getPolynomialRing().getVar(i + 1));
		}
		return AffineMorphism.fromPolynomials(domain, this, polynomials);
	}

	@Override
	public boolean isReduced() {
		return coordinateRing.isReduced();
	}

	@Override
	public boolean isIrreducible() {
		return coordinateRing.isIrreducible();
	}

	@Override
	public boolean isIntegral() {
		return coordinateRing.isIntegral();
	}

	public static class AffineSimplification<T extends Element<T>> {
		private AffineMorphism<T> simplification;
		private AffineMorphism<T> inverse;

		private AffineSimplification(AffineMorphism<T> simplification, AffineMorphism<T> inverse) {
			this.simplification = simplification;
			this.inverse = inverse;
		}

		public AffineMorphism<T> getSimplification() {
			return simplification;
		}

		public AffineMorphism<T> getInverse() {
			return inverse;
		}
	}

	public AffineSimplification<T> simplify() {
		if (simplification == null) {
			if (coordinateRing.getIdeal().generators().isEmpty()) {
				simplification = new AffineSimplification<>(identityMorphism(), identityMorphism());
				return simplification;
			}
			if (coordinateRing.getIdeal().contains(coordinateRing.getPolynomialRing().one())) {
				PolynomialRing<T> simplified = AbstractPolynomialRing.getPolynomialRing(field, 1, Monomial.GREVLEX);
				AffineScheme<T> simplifiedScheme = new AffineScheme<>(field, simplified.getUnitIdeal().divideOut());
				List<CoordinateRingElement<T>> map = new ArrayList<>();
				for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
					map.add(simplifiedScheme.getCoordinateRing().zero());
				}
				simplification = new AffineSimplification<>(
						new AffineMorphism<>(this, simplifiedScheme, Collections.singletonList(coordinateRing.zero())),
						new AffineMorphism<>(simplifiedScheme, this, map));
				return simplification;
			}
			PolynomialRing<T> polynomialRing = coordinateRing.getPolynomialRing();
			for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
				PolynomialRing<T> eliminateVariable = AbstractPolynomialRing.getPolynomialRing(field,
						new Monomial.EliminateVariableOrder(i + 1), polynomialRing.getVariableNames());
				PolynomialIdeal<T> eliminateVariableIdeal = eliminateVariable.getEmbedding(coordinateRing.getIdeal());
				Polynomial<T> elimination = eliminateVariableIdeal.generators().get(0);
				if (elimination.leadingMonomial().degree() == 1 && elimination.leadingMonomial().exponents()[i] == 1) {
					PolynomialRing<T> replacedPolynomialRing = polynomialRing.eliminateVariable(i + 1);
					int[] map = new int[polynomialRing.numberOfVariables()];
					for (int j = 0; j < polynomialRing.numberOfVariables(); j++) {
						if (j < i) {
							map[j] = j;
						} else if (j == i) {
							map[j] = -1;
						} else {
							map[j] = j - 1;
						}
					}
					Polynomial<T> replacement = replacedPolynomialRing.getEmbedding(
							eliminateVariable.subtract(eliminateVariable.getVar(i + 1), elimination), map);
					List<Polynomial<T>> mapping = new ArrayList<>();
					List<Polynomial<T>> substitutes = new ArrayList<>();
					for (int j = 0; j < polynomialRing.numberOfVariables(); j++) {
						if (i == j) {
							substitutes.add(replacement);
						} else {
							substitutes.add(replacedPolynomialRing.getEmbedding(polynomialRing.getVar(j + 1), map));
							mapping.add(polynomialRing.getVar(j + 1));
						}
					}
					List<Polynomial<T>> generators = new ArrayList<>();
					for (int j = 1; j < eliminateVariableIdeal.generators().size(); j++) {
						Polynomial<T> generator = eliminateVariableIdeal.generators().get(j);
						generators.add(replacedPolynomialRing.substitute(generator, substitutes));
					}
					AffineScheme<T> simplifiedScheme = new AffineScheme<>(field,
							replacedPolynomialRing.getIdeal(generators).divideOut());
					AffineMorphism<T> simplification = AffineMorphism.fromPolynomials(this, simplifiedScheme, mapping);
					AffineMorphism<T> inverse = AffineMorphism.fromPolynomials(simplifiedScheme, this, substitutes);
					AffineSimplification<T> furtherSimplification = simplifiedScheme.simplify();
					this.simplification = new AffineSimplification<>(
							AffineMorphism.composition(simplification, furtherSimplification.getSimplification()),
							AffineMorphism.composition(furtherSimplification.getInverse(), inverse));
					return this.simplification;
				}
			}
			simplification = new AffineSimplification<>(identityMorphism(), identityMorphism());
		}
		return simplification;
	}

	@Override
	public List<AffinePoint<T>> singularPoints() {
		List<AffinePoint<T>> result = new ArrayList<>();
		Optional<AffineMorphism<T>> singularLocus = singularLocus();
		if (singularLocus.isEmpty()) {
			return Collections.emptyList();
		}
		if (singularLocus.get().getDomain().dimension() == 0) {
			for (AffinePoint<T> point : singularLocus.get().getDomain()) {
				result.add(point);
			}
		}
//		for (AffineMorphism<T> singular : singularLocus.get().getDomain().irreducibleComponents()) {
//			if (singular.getDomain().dimension() == 0) {
//				for (AffinePoint<T> point : singular.getDomain()) {
//					result.add(point);
//				}
//			}
//		}
		return result;
	}

	@Override
	public Optional<AffineMorphism<T>> singularLocus() {
		if (coordinateRing.getIdeal().generators().size() == 0) {
			return Optional.empty();
		}
		PolynomialRing<T> r = coordinateRing.getPolynomialRing();
		DifferentialForms<T> differentialForms = r.differentialForms();
		List<Vector<Polynomial<T>>> cotangentRelations = new ArrayList<>();
		for (Polynomial<T> generator : coordinateRing.getIdeal().generators()) {
			cotangentRelations.add(differentialForms.asGradedVector(differentialForms.derivative(generator), 1));
		}
		Matrix<Polynomial<T>> cotangentMatrix = Matrix.fromColumns(cotangentRelations);
		int rank = r.numberOfVariables() - dimension();
		List<Polynomial<T>> singularLocus = new ArrayList<>();
		singularLocus.addAll(cotangentMatrix.getModule(r).minors(cotangentMatrix, rank));
		singularLocus.addAll(coordinateRing.getIdeal().generators());
		PolynomialIdeal<T> ideal = r.getIdeal(singularLocus);
		if (ideal.contains(r.one())) {
			return Optional.empty();
		}
		AffineScheme<T> domain = new AffineScheme<>(field, ideal.divideOut());
		List<Polynomial<T>> map = new ArrayList<>();
		for (int i = 0; i < r.numberOfVariables(); i++) {
			map.add(r.getVar(i + 1));
		}
		return Optional.of(AffineMorphism.fromPolynomials(domain, this, map));
	}

	public static class AffineNormalization<T extends Element<T>> {
		private AffineScheme<T> normalization;
		private AffineMorphism<T> normalizationMorphism;
		private GeneralRationalFunction<T, AffinePoint<T>, AffinePoint<T>> inverse;
		private SimpleFunctionField<T> functionField;
		private MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> functionFieldEmbedding;
		private MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> normalizationFunctionFieldEmbedding;
		private MathMap<SimpleRationalFunction<T>, TExt<T>> functionFieldEmbeddingInverse;
		private MathMap<SimpleRationalFunction<T>, TExt<T>> normalizationFunctionFieldEmbeddingInverse;

		private AffineNormalization(AffineScheme<T> normalization, AffineMorphism<T> normalizationMorphism,
				GeneralRationalFunction<T, AffinePoint<T>, AffinePoint<T>> inverse,
				SimpleFunctionField<T> functionField,
				MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> functionFieldEmbedding,
				MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> normalizationFunctionFieldEmbedding,
				MathMap<SimpleRationalFunction<T>, TExt<T>> functionFieldEmbeddingInverse,
				MathMap<SimpleRationalFunction<T>, TExt<T>> normalizationFunctionFieldEmbeddingInverse) {
			this.normalization = normalization;
			this.normalizationMorphism = normalizationMorphism;
			this.inverse = inverse;
			this.functionField = functionField;
			this.functionFieldEmbedding = functionFieldEmbedding;
			this.normalizationFunctionFieldEmbedding = normalizationFunctionFieldEmbedding;
			this.functionFieldEmbeddingInverse = functionFieldEmbeddingInverse;
			this.normalizationFunctionFieldEmbeddingInverse = normalizationFunctionFieldEmbeddingInverse;
		}

		public AffineScheme<T> getNormalization() {
			return normalization;
		}

		public AffineMorphism<T> getNormalizationMorphism() {
			return normalizationMorphism;
		}

		public GeneralRationalFunction<T, AffinePoint<T>, AffinePoint<T>> getInverse() {
			return inverse;
		}

		public SimpleFunctionField<T> getFunctionField() {
			return functionField;
		}

		public MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> getFunctionFieldEmbedding() {
			return functionFieldEmbedding;
		}

		public MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> getNormalizationFunctionFieldEmbedding() {
			return normalizationFunctionFieldEmbedding;
		}

		public MathMap<SimpleRationalFunction<T>, TExt<T>> getFunctionFieldEmbeddingInverse() {
			return functionFieldEmbeddingInverse;
		}

		public MathMap<SimpleRationalFunction<T>, TExt<T>> getNormalizationFunctionFieldEmbeddingInverse() {
			return normalizationFunctionFieldEmbeddingInverse;
		}
	}

	public AffineNormalization<T> normalization() {
//		if (!isIntegral()) {
//			throw new ArithmeticException("Not implemented!");// TODO
//		}
		if (dimension() != 1) {
			throw new ArithmeticException("Not implemented!");// TODO
		}
		return normalization(field.getExtension(field.getUnivariatePolynomialRing().getVar()));
	}

	private <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> AffineNormalization<T> normalization(
			Extension<T, B, E, FE> trivialExtension) {
		SimpleFunctionFieldFromCoordinateRing<T> functionField = SimpleFunctionField.fromCoordinateRing(coordinateRing);
		TranscendentalFieldExtension<T> transcendental = functionField.getSimpleFunctionField().getBaseField();
		Polynomial<T> denominator = transcendental.polynomialRing().one();
		for (int i = 0; i < functionField.getSimpleFunctionField().degree(); i++) {
			denominator = transcendental.polynomialRing().lcm(functionField.getSimpleFunctionField().minimalPolynomial()
					.univariateCoefficient(i).getDenominator(), denominator);
		}
		UnivariatePolynomial<TExt<T>> subs = transcendental.getUnivariatePolynomialRing()
				.getEmbedding(transcendental.getElement(transcendental.polynomialRing().one(), denominator), 1);
		UnivariatePolynomial<TExt<T>> inverseSubs = transcendental.getUnivariatePolynomialRing()
				.getEmbedding(transcendental.getEmbedding(denominator), 1);
		UnivariatePolynomial<Polynomial<T>> minimalPolynomial = transcendental.polynomialRing()
				.getUnivariatePolynomialRing().getEmbedding(
						transcendental.getUnivariatePolynomialRing()
								.normalize(transcendental.getUnivariatePolynomialRing().substitute(
										functionField.getSimpleFunctionField().minimalPolynomial(),
										Collections.singletonList(subs))),
						new FunctionMathMap<>((TExt<T> t) -> t.asInteger()));
		SimpleFunctionField<T> scaledFunctionField = new SimpleFunctionField<>(
				transcendental.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial,
						new FunctionMathMap<>((Polynomial<T> t) -> transcendental.getEmbedding(t))),
				transcendental, functionField.getSimpleFunctionField().getVariableName());
//		PolynomialRing<T> twoVariables = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.GREVLEX);
//		Polynomial<T> flattenedMinimalPolynomial = twoVariables.fromUnivariatePolynomial(minimalPolynomial, 1);
//		UnivariatePolynomial<Polynomial<T>> switchedMinimalPolynomial = twoVariables
//				.asUnivariatePolynomial(flattenedMinimalPolynomial, 2);
//		UnivariatePolynomial<TExt<T>> switchedTranscendentalMinimalPolynomial = transcendental
//				.getUnivariatePolynomialRing().getEmbedding(switchedMinimalPolynomial,
//						new FunctionMathMap<>((Polynomial<T> t) -> transcendental.getEmbedding(t)));
//		SimpleFunctionField<T> switchedFunctionField = new SimpleFunctionField<>(
//				switchedTranscendentalMinimalPolynomial, transcendental, "X");
//		MathMap<Polynomial<T>, SimpleRationalFunction<T>> switchedEmbedding = new FunctionMathMap<>(
//				(Polynomial<T> t) -> {
//					UnivariatePolynomial<Polynomial<T>> asUnivariate = twoVariables.asUnivariatePolynomial(t, 1);
//					UnivariatePolynomial<TExt<T>> asTranscendental = transcendental.getUnivariatePolynomialRing()
//							.getEmbedding(asUnivariate,
//									new FunctionMathMap<>((Polynomial<T> s) -> transcendental.getEmbedding(s)));
//					return switchedFunctionField.fromPolynomial(asTranscendental);
//				});
//		Polynomial<T> discriminant = transcendental.polynomialRing().getUnivariatePolynomialRing()
//				.discriminant(minimalPolynomial);
//		FactorizationResult<Polynomial<T>, Polynomial<T>> factors = transcendental.polynomialRing()
//				.uniqueFactorization(discriminant);
		FieldExtensionUnivariatePolynomialRing<B, E, FE> globalRing = new FieldExtensionUnivariatePolynomialRing<B, E, FE>(
				trivialExtension.extension());
		UnivariatePolynomial<Polynomial<E>> globalMinimalPolynomial = trivialExtension.extension()
				.getUnivariatePolynomialRing().getUnivariatePolynomialRing()
				.getEmbedding(minimalPolynomial, new FunctionMathMap<>((Polynomial<T> t) -> trivialExtension.extension()
						.getUnivariatePolynomialRing().getEmbedding(t, trivialExtension.embeddingMap())));
		List<UnivariatePolynomial<TExt<E>>> globalIntegralBasis = globalRing
				.extensionIntegralBasis(globalMinimalPolynomial);// new ArrayList<>();
		List<UnivariatePolynomial<TExt<T>>> integralBasis = new ArrayList<>();
		for (UnivariatePolynomial<TExt<E>> b : globalIntegralBasis) {
			integralBasis.add(transcendental.getUnivariatePolynomialRing().getEmbedding(b, new FunctionMathMap<>(
					(TExt<E> t) -> transcendental.getEmbedding(t, trivialExtension.retractionMap()))));
		}
		PolynomialRing<T> normalizationPolynomialRing = AbstractPolynomialRing.getPolynomialRing(field,
				transcendental.transcendenceDegree() + integralBasis.size()
						+ coordinateRing.getPolynomialRing().numberOfVariables(),
				Monomial.GREVLEX);
		List<Polynomial<T>> normalizationEquations = new ArrayList<>();
		List<SimpleRationalFunction<T>> asRationalFunctions = new ArrayList<>();
		List<Vector<TExt<T>>> asVectors = new ArrayList<>();
		for (UnivariatePolynomial<TExt<T>> generator : integralBasis) {
			SimpleRationalFunction<T> g = scaledFunctionField.fromPolynomial(generator);
			asRationalFunctions.add(g);
			asVectors.add(scaledFunctionField.asVector(g));
		}
		MatrixAlgebra<TExt<T>> matrixAlgebra = scaledFunctionField.matrixAlgebra();
		Matrix<TExt<T>> fromIntegralBasis = Matrix.fromColumns(asVectors);
		Matrix<TExt<T>> toIntegralBasis = matrixAlgebra.inverse(fromIntegralBasis);
		for (int i = 0; i < integralBasis.size(); i++) {
			SimpleRationalFunction<T> generator1 = asRationalFunctions.get(i);
			for (int j = 0; j <= i; j++) {
				SimpleRationalFunction<T> generator2 = asRationalFunctions.get(j);
				SimpleRationalFunction<T> product = scaledFunctionField.multiply(generator1, generator2);
				Vector<TExt<T>> productVector = matrixAlgebra.multiply(toIntegralBasis,
						scaledFunctionField.asVector(product));
				Polynomial<T> generator = normalizationPolynomialRing.multiply(-1,
						normalizationPolynomialRing.getVar(transcendental.transcendenceDegree() + i + 1),
						normalizationPolynomialRing.getVar(transcendental.transcendenceDegree() + j + 1));
				for (int k = 0; k < minimalPolynomial.degree(); k++) {
					Polynomial<T> coeff = normalizationPolynomialRing
							.getEmbedding(productVector.get(k + 1).asInteger());
					generator = normalizationPolynomialRing.add(
							normalizationPolynomialRing.multiply(coeff,
									normalizationPolynomialRing.getVar(transcendental.transcendenceDegree() + k + 1)),
							generator);
				}
				normalizationEquations.add(generator);
			}
		}
		Vector<TExt<T>> oneVector = matrixAlgebra.multiply(toIntegralBasis,
				scaledFunctionField.asVector(scaledFunctionField.one()));
		Polynomial<T> oneGenerator = normalizationPolynomialRing.negative(normalizationPolynomialRing.one());
		for (int k = 0; k < minimalPolynomial.degree(); k++) {
			Polynomial<T> coeff = normalizationPolynomialRing.getEmbedding(oneVector.get(k + 1).asInteger());
			oneGenerator = normalizationPolynomialRing.add(
					normalizationPolynomialRing.multiply(coeff,
							normalizationPolynomialRing.getVar(transcendental.transcendenceDegree() + k + 1)),
					oneGenerator);
		}
		normalizationEquations.add(oneGenerator);

		// Extension<E, B, E, FE> extension = trivialExtension.extension()
//				.getExtension(trivialExtension.extension().getUnivariatePolynomialRing().getVar());
//		for (Polynomial<T> prime : factors.primeFactors()) {
//			List<TwoGeneratorIdeal<Polynomial<E>, TExt<E>, E, B, E, FE>> ideals = globalRing.idealsOverGenerators(
//					globalMinimalPolynomial, globalRing.getEmbedding(prime, trivialExtension.embeddingMap()),
//					extension);
//			for (TwoGeneratorIdeal<Polynomial<E>, TExt<E>, E, B, E, FE> ideal : ideals) {
//				UnivariatePolynomial<TExt<T>> uniformizer = transcendental.getUnivariatePolynomialRing().getEmbedding(
//						ideal.getUniformizer(),
//						new FunctionMathMap<>((TExt<E> t) -> transcendental.getElement(
//								transcendental.polynomialRing().getEmbedding(t.getNumerator(),
//										trivialExtension.retractionMap()),
//								transcendental.polynomialRing().getEmbedding(t.getDenominator(),
//										trivialExtension.retractionMap()))));
//				TExt<T> intGenerator = transcendental.getEmbedding(transcendental.polynomialRing()
//						.getEmbedding(ideal.getIntGenerator(), trivialExtension.retractionMap()));
//				integralBasis.add(transcendental.getUnivariatePolynomialRing().divideScalar(transcendental
//						.getUnivariatePolynomialRing().power(uniformizer, ideal.getType().ramificationIndex()),
//						intGenerator));
//			}
//			LocalizedUnivariatePolynomialRing<B, E, FE> localized = new LocalizedUnivariatePolynomialRing<>(
//					trivialExtension.extension(), transcendental.polynomialRing().getVariableNames()[0],
//					trivialExtension.extension().getUnivariatePolynomialRing().getEmbedding(prime,
//							trivialExtension.embeddingMap()));
//			UnivariatePolynomial<TExt<E>> localizedMinimalPolynomial = localized.getUnivariatePolynomialRing()
//					.getEmbedding(minimalPolynomial, new FunctionMathMap<>(
//							(Polynomial<T> t) -> localized.getEmbedding(trivialExtension.extension()
//									.getUnivariatePolynomialRing().getEmbedding(t, trivialExtension.embeddingMap()))));
//			List<UnivariatePolynomial<TExt<E>>> basis = localized.ringOfIntegers().triagonalizeIntegralBasis(
//					localizedMinimalPolynomial,
//					localized.ringOfIntegers().integralBasis(localizedMinimalPolynomial,
//							localized.ringOfIntegers().theMontesAlgorithm(localizedMinimalPolynomial,
//									trivialExtension.extension().getExtension(
//											trivialExtension.extension().getUnivariatePolynomialRing().getVar())),
//							false));
//			for (int i = 0; i < functionField.getSimpleFunctionField().degree(); i++) {
//				if (localized.ringOfIntegers().valuationOfUnivariatePolynomial(basis.get(i))
//						.compareTo(Value.ZERO) < 0) {
//					integralBasis.add(transcendental.getUnivariatePolynomialRing().getEmbedding(basis.get(i),
//							new FunctionMathMap<>((TExt<E> t) -> transcendental.getElement(
//									transcendental.polynomialRing().getEmbedding(t.getNumerator(),
//											trivialExtension.retractionMap()),
//									transcendental.polynomialRing().getEmbedding(t.getDenominator(),
//											trivialExtension.retractionMap())))));
//				}
//			}
//			List<List<UnivariatePolynomial<TExt<E>>>> basis = localized.ringOfIntegers().integralBasis(
//					localizedMinimalPolynomial,
//					localized.ringOfIntegers().theMontesAlgorithm(localizedMinimalPolynomial,
//							trivialExtension.extension()
//									.getExtension(trivialExtension.extension().getUnivariatePolynomialRing().getVar())),
//					true);
//			for (List<UnivariatePolynomial<TExt<E>>> basisList : basis) {
//				for (UnivariatePolynomial<TExt<E>> integral : basisList) {
//					if (localized.ringOfIntegers().valuationOfUnivariatePolynomial(integral)
//							.compareTo(Value.ZERO) < 0) {
//						integralBasis.add(transcendental.getUnivariatePolynomialRing().getEmbedding(integral,
//								new FunctionMathMap<>((TExt<E> t) -> transcendental.getElement(
//										transcendental.polynomialRing().getEmbedding(t.getNumerator(),
//												trivialExtension.retractionMap()),
//										transcendental.polynomialRing().getEmbedding(t.getDenominator(),
//												trivialExtension.retractionMap())))));
//					}
//				}
//			}
//		}

//		PolynomialRing<T> normalizationPolynomialRing = AbstractPolynomialRing.getPolynomialRing(field,
//				transcendental.transcendenceDegree() + integralBasis.size() + 1
//						+ coordinateRing.getPolynomialRing().numberOfVariables(),
//				Monomial.GREVLEX);
//		List<Polynomial<T>> normalizationEquations = new ArrayList<>();
//		Polynomial<T> embeddedMinimalPolynomial = normalizationPolynomialRing.zero();
//		for (int i = 0; i <= minimalPolynomial.degree(); i++) {
//			embeddedMinimalPolynomial = normalizationPolynomialRing.add(
//					normalizationPolynomialRing.multiply(
//							normalizationPolynomialRing.getEmbedding(minimalPolynomial.univariateCoefficient(i)),
//							normalizationPolynomialRing.getVarPower(transcendental.transcendenceDegree() + 1, i)),
//					embeddedMinimalPolynomial);
//		}
//		normalizationEquations.add(embeddedMinimalPolynomial);
//		CoordinateRing<T> normalizationCoordinateRing = normalizationPolynomialRing.getIdeal(normalizationEquations)
//				.divideOut();
//		for (int j = 0; j < integralBasis.size(); j++) {
//			normalizationCoordinateRing = normalizationPolynomialRing.getIdeal(normalizationEquations).divideOut();
//			normalizationEquations.clear();
//			normalizationEquations.addAll(normalizationCoordinateRing.getIdeal().generators());
//			UnivariatePolynomial<TExt<T>> integral = integralBasis.get(j);
//			normalizationEquations.addAll(integerMinimalPolynomials(normalizationCoordinateRing, integral,
//					transcendental.transcendenceDegree() + j + 2, functionField.getSimpleFunctionField()));
//			SimpleRationalFunction<T> morphism = functionField.getSimpleFunctionField().fromPolynomial(integral);
//			UnivariatePolynomial<TExt<T>> integerPolynomial = functionField.getSimpleFunctionField()
//					.minimalPolynomial(morphism);
//			embeddedMinimalPolynomial = normalizationPolynomialRing.zero();
//			for (int i = 0; i <= integerPolynomial.degree(); i++) {
//				embeddedMinimalPolynomial = normalizationPolynomialRing.add(normalizationPolynomialRing.multiply(
//						normalizationPolynomialRing
//								.getEmbedding(integerPolynomial.univariateCoefficient(i).asInteger()),
//						normalizationPolynomialRing.getVarPower(transcendental.transcendenceDegree() + j + 2, i)),
//						embeddedMinimalPolynomial);
//			}
//			normalizationEquations.add(embeddedMinimalPolynomial);
//			Polynomial<T> basisDenominator = transcendental.polynomialRing().one();
//			for (int i = 0; i < minimalPolynomial.degree(); i++) {
//				basisDenominator = transcendental.polynomialRing()
//						.lcm(integral.univariateCoefficient(i).getDenominator(), basisDenominator);
//			}
//			Polynomial<T> basisPolynomial = normalizationPolynomialRing.multiply(-1,
//					normalizationPolynomialRing.getEmbedding(basisDenominator),
//					normalizationPolynomialRing.getVar(transcendental.transcendenceDegree() + j + 2));
//			Polynomial<T> numerator = twoVariables.zero();
//			for (int i = 0; i < minimalPolynomial.degree(); i++) {
//				Polynomial<T> numeratorCoeff = twoVariables.getEmbedding(transcendental
//						.multiply(transcendental.getEmbedding(basisDenominator), integral.univariateCoefficient(i))
//						.asInteger());
//				numeratorCoeff = twoVariables.multiply(
//						twoVariables.getVarPower(transcendental.transcendenceDegree() + 1, i), numeratorCoeff);
//				numerator = twoVariables.add(numerator, numeratorCoeff);
//			}
//			basisPolynomial = normalizationPolynomialRing.add(basisPolynomial,
//					normalizationPolynomialRing.getEmbedding(numerator));
//			normalizationEquations.add(basisPolynomial);
//
//			SimpleRationalFunction<T> asInverseRationalFunction = switchedFunctionField.divide(
//					switchedEmbedding.evaluate(numerator),
//					switchedEmbedding.evaluate(twoVariables.getEmbedding(basisDenominator)));
//			int[] invertVars = new int[] { 1, 0 };
//			UnivariatePolynomial<TExt<T>> inverse = asInverseRationalFunction.asPolynomial();
//			Polynomial<T> inverseDenominator = transcendental.polynomialRing().one();
//			for (int i = 0; i < switchedFunctionField.degree(); i++) {
//				inverseDenominator = transcendental.polynomialRing()
//						.lcm(inverse.univariateCoefficient(i).getDenominator(), inverseDenominator);
//			}
//			Polynomial<T> inversePolynomial = normalizationPolynomialRing.multiply(-1,
//					normalizationPolynomialRing.getEmbedding(inverseDenominator, invertVars),
//					normalizationPolynomialRing.getVar(transcendental.transcendenceDegree() + j + 2));
//			for (int i = 0; i < switchedFunctionField.degree(); i++) {
//				numerator = normalizationPolynomialRing.getEmbedding(transcendental
//						.multiply(transcendental.getEmbedding(inverseDenominator), inverse.univariateCoefficient(i))
//						.asInteger(), invertVars);
//				numerator = normalizationPolynomialRing.multiply(
//						normalizationPolynomialRing.getVarPower(transcendental.transcendenceDegree(), i), numerator);
//				inversePolynomial = normalizationPolynomialRing.add(inversePolynomial, numerator);
//			}
//			normalizationEquations.add(inversePolynomial);
//			integerPolynomial = switchedFunctionField.minimalPolynomial(asInverseRationalFunction);
//			embeddedMinimalPolynomial = normalizationPolynomialRing.zero();
//			for (int i = 0; i <= integerPolynomial.degree(); i++) {
//				embeddedMinimalPolynomial = normalizationPolynomialRing.add(normalizationPolynomialRing.multiply(
//						normalizationPolynomialRing.getEmbedding(integerPolynomial.univariateCoefficient(i).asInteger(),
//								invertVars),
//						normalizationPolynomialRing.getVarPower(transcendental.transcendenceDegree() + j + 2, i)),
//						embeddedMinimalPolynomial);
//			}
//			normalizationEquations.add(embeddedMinimalPolynomial);
//		}
		Vector<TExt<T>> generatorVector = matrixAlgebra.multiply(toIntegralBasis,
				scaledFunctionField.asVector(scaledFunctionField.alpha()));
		Polynomial<T> generator = normalizationPolynomialRing.zero();
		for (int k = 0; k < minimalPolynomial.degree(); k++) {
			Polynomial<T> coeff = normalizationPolynomialRing.getEmbedding(generatorVector.get(k + 1).asInteger());
			generator = normalizationPolynomialRing.add(
					normalizationPolynomialRing.multiply(coeff,
							normalizationPolynomialRing.getVar(transcendental.transcendenceDegree() + k + 1)),
					generator);
		}
		List<Polynomial<T>> mapping = new ArrayList<>();
		for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
			SimpleRationalFunction<T> asRationalFunction = functionField.getIsomorphism()
					.evaluate(coordinateRing.getVar(i + 1));
			UnivariatePolynomial<TExt<T>> scaled = transcendental.getUnivariatePolynomialRing()
					.substitute(asRationalFunction.asPolynomial(), Collections.singletonList(subs));
			Polynomial<T> basisDenominator = transcendental.polynomialRing().one();
			for (int j = 0; j < minimalPolynomial.degree(); j++) {
				basisDenominator = transcendental.polynomialRing().lcm(scaled.univariateCoefficient(j).getDenominator(),
						basisDenominator);
			}
			Polynomial<T> basisPolynomial = normalizationPolynomialRing.multiply(-1,
					normalizationPolynomialRing.getEmbedding(basisDenominator), normalizationPolynomialRing
							.getVar(transcendental.transcendenceDegree() + integralBasis.size() + i + 1));
			for (int j = 0; j < minimalPolynomial.degree(); j++) {
				Polynomial<T> numerator = normalizationPolynomialRing.getEmbedding(transcendental
						.multiply(transcendental.getEmbedding(basisDenominator), scaled.univariateCoefficient(j))
						.asInteger());
				numerator = normalizationPolynomialRing.multiply(normalizationPolynomialRing.power(generator, j),
						numerator);
				basisPolynomial = normalizationPolynomialRing.add(basisPolynomial, numerator);
			}
			normalizationEquations.add(basisPolynomial);
			mapping.add(normalizationPolynomialRing
					.getVar(transcendental.transcendenceDegree() + integralBasis.size() + i + 1));
		}
		for (Polynomial<T> rangeRelation : coordinateRing.getIdeal().generators()) {
			normalizationEquations.add(normalizationPolynomialRing.substitute(rangeRelation, mapping));
		}
		AffineScheme<T> unsimplifiedNormalization = new AffineScheme<>(field,
				normalizationPolynomialRing.getIdeal(normalizationEquations).divideOut());
		AffineMorphism<T> unsimplifiedNormalizationMorphism = AffineMorphism.fromPolynomials(unsimplifiedNormalization,
				this, mapping);
		AffineSimplification<T> simplification = unsimplifiedNormalization.simplify();
		AffineScheme<T> normalization = simplification.getSimplification().getRange();
		AffineMorphism<T> normalizationMorphism = AffineMorphism.composition(simplification.getInverse(),
				unsimplifiedNormalizationMorphism);
		MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> coordinateRingToFunctionFieldMap = new FunctionMathMap<>(
				(CoordinateRingElement<T> t) -> scaledFunctionField
						.fromPolynomial(transcendental.getUnivariatePolynomialRing().substitute(
								functionField.getIsomorphism().evaluate(t).asPolynomial(),
								Collections.singletonList(subs))));
		MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> normalizationCoordinateRingToFunctionFieldMap = new FunctionMathMap<>(
				(CoordinateRingElement<T> t) -> {
					CoordinateRingElement<T> s = simplification.getSimplification().inducedMap(t);
					UnivariatePolynomialRing<TExt<T>> univariate = transcendental.getUnivariatePolynomialRing();
					List<Polynomial<TExt<T>>> substitutes = new ArrayList<>();
					substitutes.add(univariate.getEmbedding(transcendental.getVar(1)));
					for (Polynomial<TExt<T>> integral : integralBasis) {
						substitutes.add(integral);
					}
					for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
						substitutes.add(
								coordinateRingToFunctionFieldMap.evaluate(coordinateRing.getVar(i + 1)).asPolynomial());
					}
					UnivariatePolynomial<TExt<T>> asPolynomial = univariate.zero();
					for (Monomial m : s.getElement().monomials()) {
						UnivariatePolynomial<TExt<T>> coeff = univariate
								.getEmbedding(transcendental.getEmbedding(s.getElement().coefficient(m)));
						int[] exp = m.exponents();
						for (int i = 0; i < exp.length; i++) {
							coeff = univariate.multiply(coeff, univariate.power(substitutes.get(i), exp[i]));
						}
						asPolynomial = univariate.add(coeff, asPolynomial);
					}
					return scaledFunctionField.fromPolynomial(asPolynomial);
				});
		TranscendentalFieldExtension<T> normalizationTranscendental = new TranscendentalFieldExtension<T>(field,
				normalization.getCoordinateRing().getPolynomialRing());
		MathMap<SimpleRationalFunction<T>, TExt<T>> functionFieldToNormalizationCoordinateRingMap = new FunctionMathMap<>(
				(SimpleRationalFunction<T> t) -> {
					Vector<TExt<T>> asVector = scaledFunctionField.asVector(t);
					Vector<TExt<T>> inIntegralBasis = matrixAlgebra.multiply(toIntegralBasis, asVector);
					Polynomial<T> denominatorPolynomial = transcendental.polynomialRing().one();
					for (TExt<T> coeff : inIntegralBasis.asList()) {
						denominatorPolynomial = transcendental.polynomialRing().lcm(coeff.getDenominator(),
								denominatorPolynomial);
					}
					Polynomial<T> embeddedDenominator = normalizationPolynomialRing.getEmbedding(denominatorPolynomial);
					Polynomial<T> numeratorPolynomial = normalizationPolynomialRing.zero();
					for (int i = 0; i < scaledFunctionField.degree(); i++) {
						TExt<T> coeff = inIntegralBasis.get(i + 1);
						Polynomial<T> numerator = normalizationPolynomialRing.getEmbedding(transcendental
								.multiply(transcendental.getEmbedding(denominatorPolynomial), coeff).asInteger());
						numeratorPolynomial = normalizationPolynomialRing.add(
								normalizationPolynomialRing.multiply(numerator,
										normalizationPolynomialRing
												.getVar(transcendental.transcendenceDegree() + i + 1)),
								numeratorPolynomial);
					}
					return normalizationTranscendental.getElement(
							simplification.getInverse().inducedMap(numeratorPolynomial),
							simplification.getInverse().inducedMap(embeddedDenominator));
				});
		MathMap<SimpleRationalFunction<T>, TExt<T>> functionFieldToCoordinateRingMap = new FunctionMathMap<>(
				(SimpleRationalFunction<T> t) -> {
					UnivariatePolynomial<TExt<T>> asPolynomial = t.asPolynomial();
					UnivariatePolynomial<TExt<T>> backSubstituted = transcendental.getUnivariatePolynomialRing()
							.substitute(asPolynomial, Collections.singletonList(inverseSubs));
					return functionField.getInverseIsomorphism()
							.evaluate(functionField.getSimpleFunctionField().fromPolynomial(backSubstituted));
				});
		List<TExt<T>> toNormalization = new ArrayList<>();
		Polynomial<T> integralDenominator = coordinateRing.getPolynomialRing().one();
		for (int i = 0; i < normalization.getCoordinateRing().getPolynomialRing().numberOfVariables(); i++) {
			CoordinateRingElement<T> var = normalization.getCoordinateRing().getVar(i + 1);
			SimpleRationalFunction<T> asRationalFunction = normalizationCoordinateRingToFunctionFieldMap.evaluate(var);
			TExt<T> asFraction = functionFieldToCoordinateRingMap.evaluate(asRationalFunction);
			toNormalization.add(asFraction);
			integralDenominator = coordinateRing.getPolynomialRing().lcm(asFraction.getDenominator(),
					integralDenominator);
		}
		OpenImmersionResult<T> restriction = AffineMorphism.getOpenImmersion(this, integralDenominator);
		integralDenominator = restriction.getPolynomial().getElement();
		TranscendentalFieldExtension<T> baseTranscendental = new TranscendentalFieldExtension<>(field,
				coordinateRing.getPolynomialRing());
		PolynomialRing<T> restrictionPolynomialRing = restriction.getOpenImmersion().getDomain().getCoordinateRing()
				.getPolynomialRing();
		List<Polynomial<T>> toNormalizationMorphism = new ArrayList<>();
		for (int i = 0; i < toNormalization.size(); i++) {
			int power = 0;
			TExt<T> image = toNormalization.get(i);
			while (!coordinateRing.getPolynomialRing().isUnit(image.getDenominator())) {
				power++;
				image = baseTranscendental.multiply(baseTranscendental.getEmbedding(integralDenominator), image);
			}
			toNormalizationMorphism.add(restrictionPolynomialRing.multiply(
					restrictionPolynomialRing.power(restriction.getInverse().getElement(), power),
					restriction.getOpenImmersion().inducedMap(image.asInteger())));
		}
		AffineMorphism<T> toNormalziationMap = AffineMorphism
				.fromPolynomials(restriction.getOpenImmersion().getDomain(), normalization, toNormalizationMorphism);
		GeneralRationalFunction<T, AffinePoint<T>, AffinePoint<T>> inverse = new GeneralRationalFunction<>(this,
				normalization, toNormalziationMap, 0, restriction.getOpenImmersion(), 0);
		return new AffineNormalization<T>(normalization, normalizationMorphism, inverse, scaledFunctionField,
				coordinateRingToFunctionFieldMap, normalizationCoordinateRingToFunctionFieldMap,
				functionFieldToCoordinateRingMap, functionFieldToNormalizationCoordinateRingMap);
	}

	private Set<Polynomial<T>> integerMinimalPolynomials(CoordinateRing<T> normalizationCoordinateRing,
			UnivariatePolynomial<TExt<T>> integral, int var, SimpleFunctionField<T> functionField) {
		PolynomialRing<T> normalizationPolynomialRing = normalizationCoordinateRing.getPolynomialRing();
		TranscendentalFieldExtension<T> transcendental = functionField.getBaseField();
//		PolynomialRing<T> basePolynomialRing = transcendental.polynomialRing();
//		Polynomial<T> denominator = basePolynomialRing.one();
//		for (int i = 0; i <= integral.degree(); i++) {
//			denominator = basePolynomialRing.lcm(integral.univariateCoefficient(i).getDenominator(), denominator);
//		}
//		PolynomialRing<T> added = AbstractPolynomialRing.getPolynomialRing(field,
//				normalizationPolynomialRing.numberOfVariables() + 1,
//				new Monomial.InvertedEliminateVariableOrder(normalizationPolynomialRing.numberOfVariables() + 1));
//		PolynomialIdeal<T> ideal = added.getEmbedding(normalizationCoordinateRing.getIdeal());
//		ideal = added
//				.add(ideal,
//						added.getIdeal(
//								Collections
//										.singletonList(added.subtract(
//												added.multiply(added.getEmbedding(denominator),
//														added.getVar(
//																normalizationPolynomialRing.numberOfVariables() + 1)),
//												added.one()))));
//		CoordinateRing<T> cr = ideal.divideOut();
//		CoordinateRingElement<T> numerator = cr.zero();
//		for (int i = 0; i <= integral.degree(); i++) {
//			numerator = cr.add(
//					cr.multiply(
//							cr.getEmbedding(
//									added.getEmbedding(transcendental.multiply(integral.univariateCoefficient(i),
//											transcendental.getEmbedding(denominator)).asInteger())),
//							cr.getEmbedding(added.getVarPower(transcendental.transcendenceDegree() + 1, i))),
//					numerator);
//		}
//		CoordinateRingElement<T> integerElement = cr
//				.multiply(cr.getVar(normalizationPolynomialRing.numberOfVariables() + 1), numerator);
//		List<CoordinateRingElement<T>> powers = new ArrayList<>();
//		powers.add(cr.one());
//		CoordinateRingElement<T> power = cr.one();
//		Set<Polynomial<T>> result = new TreeSet<>();
//		while (result.isEmpty() && powers.size() <= functionField.degree()) {
//			power = cr.multiply(integerElement, power);
//			powers.add(power);
//			//X^3 + (-2*Y + 2)*X^2 + -1*Y^2*X + (2*Y^5 + 2*Y^4)
//			List<Vector<CoordinateRingElement<T>>> syzygies = cr.syzygyProblem(powers);
//			for (Vector<CoordinateRingElement<T>> syzygy : syzygies) {
//				boolean valid = true;
//				for (CoordinateRingElement<T> coeff : syzygy.asList()) {
//					if (!coeff.equals(cr.zero())
//							&& coeff.getElement().degree(normalizationPolynomialRing.numberOfVariables() + 1) != 0) {
//						valid = false;
//						break;
//					}
//				}
//				if (!valid) {
//					continue;
//				}
//				Polynomial<T> syzygyPolynomial = normalizationPolynomialRing.zero();
//				for (int i = 0; i < powers.size(); i++) {
//					syzygyPolynomial = normalizationPolynomialRing.add(normalizationPolynomialRing.multiply(
//							normalizationPolynomialRing.getEmbedding(syzygy.get(i + 1).getElement()),
//							normalizationPolynomialRing.getVarPower(var, i)), syzygyPolynomial);
//				}
//				result.add(syzygyPolynomial);
//			}
//		}
//		if (!result.isEmpty()) {
//			return result;
//		}
		SimpleRationalFunction<T> morphism = functionField.fromPolynomial(integral);
		UnivariatePolynomial<TExt<T>> integerPolynomial = functionField.minimalPolynomial(morphism);
		Polynomial<T> embeddedMinimalPolynomial = normalizationPolynomialRing.zero();
		for (int i = 0; i <= integerPolynomial.degree(); i++) {
			embeddedMinimalPolynomial = normalizationPolynomialRing.add(normalizationPolynomialRing.multiply(
					normalizationPolynomialRing.getEmbedding(integerPolynomial.univariateCoefficient(i).asInteger()),
					normalizationPolynomialRing.getVarPower(var, i)), embeddedMinimalPolynomial);
		}
		return Collections.singleton(embeddedMinimalPolynomial);
	}
}
