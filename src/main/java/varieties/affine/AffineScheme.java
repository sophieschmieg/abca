package varieties.affine;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.DifferentialForms;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import varieties.AbstractScheme;
import varieties.Morphism;
import varieties.SpectrumOfField;
import varieties.SpectrumOfField.SingletonPoint;

public class AffineScheme<T extends Element<T>> extends AbstractScheme<T, AffinePoint<T>> {
	private CoordinateRing<T> coordinateRing;
	private Field<T> field;
	private AffineCover<T> cover;

	public AffineScheme(Field<T> field, CoordinateRing<T> coordinateRing) {
		if (!coordinateRing.getPolynomialRing().getRing().equals(field)) {
			throw new ArithmeticException("Incompatible parameters");
		}
		this.coordinateRing = coordinateRing;
		this.field = field;
		this.cover = new AffineCover<>(Collections.singletonList(this),
				Collections.singletonList(Collections.singletonList(identityMorphism())));
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

	public static <T extends Element<T>> AffineMorphism<T> restrictAwayFrom(AffineScheme<T> variety,
			List<CoordinateRingElement<T>> constraints) {
		int numVars = variety.coordinateRing.getPolynomialRing().numberOfVariables();
		PolynomialRing<T> polynomialRing = AbstractPolynomialRing.getPolynomialRing(variety.field,
				numVars + constraints.size(), variety.coordinateRing.getPolynomialRing().getComparator());
		List<Polynomial<T>> idealGenerators = new ArrayList<>();
		for (Polynomial<T> generator : variety.coordinateRing.getIdeal().generators()) {
			idealGenerators.add(polynomialRing.getEmbedding(generator));
		}
		for (int i = 0; i < constraints.size(); i++) {
			Polynomial<T> t = polynomialRing.getEmbedding(constraints.get(i).getElement());
			idealGenerators.add(polynomialRing.subtract(
					polynomialRing.multiply(t, polynomialRing.getVar(i + 1 + numVars)), polynomialRing.one()));
		}
		PolynomialIdeal<T> ideal = polynomialRing.getIdeal(idealGenerators);
		List<Polynomial<T>> identity = new ArrayList<>();
		for (int i = 0; i < numVars; i++) {
			identity.add(polynomialRing.getVar(i + 1));
		}
		AffineScheme<T> restricted = new AffineScheme<>(variety.getField(), ideal.divideOut());
		return AffineMorphism.fromPolynomials(restricted, variety, identity);
	}

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
		CoordinateIdeal<T> ideal = coordinateRing.getIdeal(point.asIdeal(coordinateRing.getPolynomialRing()));
		return new LocalizedCoordinateRing<>(field, coordinateRing, ideal);
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
	public AffinePoint<T> getRandomElement() {
		throw new UnsupportedOperationException();
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
		return new Morphism<>() {

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
			public RestrictionResult<T> restrict(SingletonPoint preimage) {
				AffineScheme<T> singleton = domain.getAffineCover().getCover().get(0);
				PolynomialRing<T> polynomials = singleton.coordinateRing.getPolynomialRing();
				List<Polynomial<T>> asPolynomials = new ArrayList<>();
				for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
					asPolynomials.add(polynomials.getEmbedding(p.getCoord(i + 1)));
				}
				return new RestrictionResult<>(0, 0, singleton.identityMorphism(), identityMorphism(),
						AffineMorphism.fromPolynomials(singleton, AffineScheme.this, asPolynomials));
			}

		};
	}

	@Override
	public List<AffineScheme<T>> irreducibleComponents() {
		List<AffineScheme<T>> result = new ArrayList<>();
		List<PolynomialIdeal<T>> components = coordinateRing.getPolynomialRing()
				.primaryDecomposition(coordinateRing.getIdeal()).getPrimaryIdeals();
		for (PolynomialIdeal<T> component : components) {
			result.add(new AffineScheme<>(field, component.divideOut()));
		}
		return result;
	}

	@Override
	public AffineScheme<T> reduced() {
		return new AffineScheme<>(field,
				coordinateRing.getPolynomialRing().radical(coordinateRing.getIdeal()).divideOut());
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

	@Override
	public List<AffinePoint<T>> singularPoints() {
		List<AffinePoint<T>> result = new ArrayList<>();
		Optional<AffineMorphism<T>> singularLocus = singularLocus();
		if (singularLocus.isEmpty()) {
			return Collections.emptyList();
		}
		for (AffineScheme<T> singular : singularLocus.get().getDomain().irreducibleComponents()) {
			if (singular.dimension() == 0) {
				for (AffinePoint<T> point : singular) {
					result.add(point);
				}
			}
		}
		return result;
	}

	@Override
	public Optional<AffineMorphism<T>> singularLocus() {
		PolynomialRing<T> r = coordinateRing.getPolynomialRing();
		DifferentialForms<T> differentialForms = new DifferentialForms<>(r);
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
}
