package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Module;
import fields.interfaces.Ring;
import fields.interfaces.Ring.ModuloIdealResult;
import fields.vectors.GenericPIDModule.Mod;

public class GenericPIDModule<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, Mod<S>>
		implements Module<T, Mod<S>> {
	public static class Mod<S extends Element<S>> extends AbstractElement<Mod<S>> implements Element<Mod<S>> {
		private S representative;

		private Mod(S representative) {
			this.representative = representative;
		}

		@Override
		public int compareTo(Mod<S> o) {
			return representative.compareTo(o.representative);
		}

		@Override
		public String toString() {
			return representative.toString();
		}
	}

	private Ring<T> base;
	private Matrix<T> syzygiesMatrix;
	private MatrixModule<T> matrixModule;
	private Module<T, S> freeModule;
	private FreeSubModule<T, S> syzygies;
	private List<Ideal<T>> diagonalIdeals;
	private List<Mod<S>> generators;
	private int embeddingDimension;

	public GenericPIDModule(Module<T, S> freeModule, FreeSubModule<T, S> syzygies) {
		this.freeModule = freeModule;
		this.base = freeModule.getRing();
		this.embeddingDimension = freeModule.getModuleGenerators().size();
		this.syzygies = syzygies;
		List<Vector<T>> vectorSyzygies = new ArrayList<>();
		for (S basisVector : syzygies.getBasis()) {
			vectorSyzygies.add(freeModule.asVector(basisVector));
		}
//		if (!vectorSyzygies.isEmpty() && freeModule.getRing() instanceof Integers) {
//			List<Vector<IntE>> integerSyzygies = new ArrayList<>();
//			MathMap<T, IntE> castMap = new MathMap<>() {
//				@Override
//				public IntE evaluate(T t) {
//					return (IntE) t;
//				}
//			};
//			for (S syzygyVector : syzygies.getBasis()) {
//				integerSyzygies.add(Vector.mapVector(castMap, freeModule.asVector(syzygyVector)));
//			}
//			RealLattice asRealLattice = RealLattice.fromIntegerLattice(Reals.r(128), embeddingDimension,
//					integerSyzygies, 1.0);
//			vectorSyzygies.clear();
//			integerSyzygies.clear();
//			List<S> syzygiesList = new ArrayList<>();
//			MathMap<IntE, T> castBackMap = new MathMap<>() {
//				@SuppressWarnings("unchecked")
//				@Override
//				public T evaluate(IntE t) {
//					return (T) t;
//				}
//			};
//			MathMap<Real, IntE> roundMap = new MathMap<>() {
//				@Override
//				public IntE evaluate(Real t) {
//					return t.round();
//				}
//			};
//			for (Vector<Real> reduced : asRealLattice.getModuleGenerators()) {
//				Vector<IntE> integerVector = Vector.mapVector(roundMap, reduced);
//				Vector<T> vectorSyzygy = Vector.mapVector(castBackMap, integerVector);
//				S syzygy = freeModule.fromVector(vectorSyzygy);
//				integerSyzygies.add(integerVector);
//				vectorSyzygies.add(vectorSyzygy);
//				syzygiesList.add(syzygy);
//			}
//			this.syzygies = new FreeSubModule<>(freeModule, syzygiesList);
//			this.rationalLattice = RationalLattice.fromIntegerLattice(embeddingDimension, integerSyzygies);
//			this.asFractionMap = new MathMap<>() {
//				@Override
//				public Fraction evaluate(T t) {
//					return Rationals.q().getEmbedding((IntE) t);
//				}
//			};
//			this.fromFractionMap = new MathMap<>() {
//				@SuppressWarnings("unchecked")
//				@Override
//				public T evaluate(Fraction t) {
//					return (T) t.asInteger();
//				}
//			};
//		}
		if (vectorSyzygies.isEmpty()) {
			vectorSyzygies.add(freeModule.asVector(freeModule.zero()));
		}
		this.syzygiesMatrix = Matrix.fromColumns(vectorSyzygies);
		this.matrixModule = new MatrixModule<>(freeModule.getRing(), embeddingDimension, vectorSyzygies.size());
		generators = new ArrayList<>();
		for (S basisVector : freeModule.getModuleGenerators()) {
			generators.add(reduce(basisVector));
		}
	}

	public static <T extends Element<T>> GenericPIDModule<T, Vector<T>> fromSyzygies(FreeModule<T> freeModule,
			List<List<T>> syzygies) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (List<T> syzygy : syzygies) {
			asVectors.add(new Vector<>(syzygy));
		}
		FreeSubModule<T, Vector<T>> subModule = new FreeSubModule<>(freeModule, asVectors);
		return new GenericPIDModule<>(freeModule, subModule);
	}

	public static <T extends Element<T>> GenericPIDModule<T, Vector<T>> fromSyzygies(Ring<T> base, int rank,
			List<List<T>> syzygies) {
		return fromSyzygies(new FreeModule<>(base, rank), syzygies);
	}

	@Override
	public String toString() {
		return freeModule + "/" + syzygies;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	public Module<T, S> getFreeModule() {
		return freeModule;
	}

	@Override
	public Mod<S> getRandomElement() {
		List<T> asVector = new ArrayList<>();
		for (int i = 0; i < embeddingDimension; i++) {
			T order = diagonalIdeals().get(i).generators().get(0);
			if (order.equals(base.zero())) {
				asVector.add(base.getRandomElement());
			} else {
				asVector.add(base.getRandomElement(order));
			}
		}
		return fromDiagonalVector(new Vector<>(asVector));
	}

	@Override
	public boolean isFinite() {
		if (base.isFinite()) {
			return true;
		}
		if (embeddingDimension != syzygies.rank()) {
			return base.isFinite();
		}
		for (Ideal<T> diagonalIdeal : diagonalIdeals()) {
			if (!base.moduloIdeal(diagonalIdeal).getQuotientRing().isFinite()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		if (!isFinite()) {
			throw new InfinityException();
		}
		BigInteger result = BigInteger.ONE;
		for (Ideal<T> diagonalIdeal : diagonalIdeals()) {
			result = result.multiply(base.moduloIdeal(diagonalIdeal).getQuotientRing().getNumberOfElements());
		}
		return result;
	}

	private <Q extends Element<Q>> Iterator<T> moduloIterator(ModuloIdealResult<T, Q> modulo) {
		return new Iterator<>() {
			private Iterator<Q> it = modulo.getQuotientRing().iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public T next() {
				return modulo.lift(it.next());
			}
		};
	}

	@Override
	public Iterator<Mod<S>> iterator() {
		List<Ideal<T>> ideals = diagonalIdeals();
		List<ModuloIdealResult<T, ?>> modulo = new ArrayList<>();
		for (Ideal<T> ideal : ideals) {
			modulo.add(base.moduloIdeal(ideal));
		}
		return new Iterator<>() {
			private List<Iterator<T>> it = null;
			private List<T> list = new ArrayList<T>();

			private void init() {
				if (this.it != null) {
					return;
				}
				this.it = new ArrayList<Iterator<T>>();
				for (int i = 0; i < embeddingDimension; i++) {
					this.it.add(moduloIterator(modulo.get(i)));
					if (i == 0) {
						this.list.add(null);
					} else {
						this.list.add(this.it.get(i).next());
					}
				}
			}

			@Override
			public boolean hasNext() {
				init();
				for (int i = 0; i < embeddingDimension; i++) {
					if (this.it.get(i).hasNext()) {
						return true;
					}
				}
				return false;
			}

			@Override
			public Mod<S> next() {
				init();
				boolean broken = false;
				for (int i = 0; i < embeddingDimension; i++) {
					if (this.it.get(i).hasNext()) {
						this.list.set(i, this.it.get(i).next());
						for (int j = 0; j < i; j++) {
							this.it.set(j, moduloIterator(modulo.get(j)));
							this.list.set(j, this.it.get(j).next());
						}
						broken = true;
						break;
					}
				}
				if (!broken) {
					throw new RuntimeException();
				}
				Mod<S> result = fromDiagonalVector(new Vector<>(list));
				return result;
			}
		};
	}

	@Override
	public Ring<T> getRing() {
		return base;
	}

	@Override
	public Mod<S> zero() {
		return reduce(freeModule.zero());
	}

	@Override
	public Mod<S> add(Mod<S> s1, Mod<S> s2) {
		return reduce(freeModule.add(s1.representative, s2.representative));
	}

	@Override
	public Mod<S> negative(Mod<S> s) {
		return reduce(freeModule.negative(s.representative));
	}

	@Override
	public Mod<S> scalarMultiply(T t, Mod<S> s) {
		return reduce(freeModule.scalarMultiply(t, s.representative));
	}

	@Override
	public boolean isFree() {
		return syzygies.rank() == 0;
	}

	@Override
	public Ideal<T> annihilator() {
		return diagonalIdeals().get(embeddingDimension - 1);
	}

	public T index() {
		T index = base.one();
		for (Ideal<T> diagonalIdeal : diagonalIdeals()) {
			index = base.multiply(diagonalIdeal.generators().get(0), index);
		}
		return index;
	}

	@Override
	public boolean isLinearIndependent(List<Mod<S>> s) {
		List<S> vectors = new ArrayList<>();
		for (Mod<S> v : s) {
			vectors.add(v.representative);
			Vector<T> vector = asDiagonalVector(v);
			for (int i = 0; i < embeddingDimension; i++) {
				T coeff = vector.get(i + 1);
				if (coeff.equals(base.zero())) {
					continue;
				}
				if (diagonalIdeals().get(i).equals(base.getZeroIdeal())) {
					continue;
				}
				return false;
			}
		}
		return freeModule.isLinearIndependent(vectors);
	}

	@Override
	public boolean isGeneratingModule(List<Mod<S>> s) {
		List<S> lifts = new ArrayList<>();
		for (Mod<S> v : s) {
			lifts.add(lift(v));
		}
		return freeModule.isGeneratingModule(lifts);
	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<Mod<S>> s) {
		List<S> lifts = new ArrayList<>();
		List<Vector<T>> result = new ArrayList<>();
		List<Ideal<T>> diagonalIdeals = diagonalIdeals();
		for (int k = 0; k < s.size(); k++) {
			Mod<S> v = s.get(k);
			lifts.add(lift(v));
			Vector<T> asDiagonalVector = asDiagonalVector(v);
			T lcm = base.one();
			for (int i = 0; i < embeddingDimension; i++) {
				T diagonalGenerator = diagonalIdeals.get(i).generators().get(0);
				T gcd = base.gcd(asDiagonalVector.get(i + 1), diagonalGenerator);
				T multiplier = base.divideChecked(diagonalGenerator, gcd);
				lcm = base.lcm(lcm, multiplier);
			}
			List<T> additionalRelation = new ArrayList<>();
			if (!lcm.equals(base.zero())) {
				for (int i = 0; i < k; i++) {
					additionalRelation.add(base.zero());
				}
				additionalRelation.add(lcm);
				for (int i = k + 1; i < s.size(); i++) {
					additionalRelation.add(base.zero());
				}
			}
			result.add(new Vector<>(additionalRelation));
		}
		result.addAll(freeModule.nonTrivialCombinations(lifts));

		return result;
	}

	@Override
	public List<Mod<S>> getModuleGenerators() {
		return generators;
	}

	@Override
	public List<Vector<T>> getSyzygies() {
		List<Vector<T>> result = new ArrayList<>();
		for (S syzygy : syzygies.getBasis()) {
			result.add(freeModule.asVector(syzygy));
		}
		return result;
	}

	@Override
	public Vector<T> asVector(Mod<S> s) {
		return freeModule.asVector(lift(s));
	}

	public List<T> diagonalRanks() {
		Matrix<T> diagonal = matrixModule.smithNormalForm(syzygiesMatrix).getDiagonalMatrix();
		List<T> result = new ArrayList<>();
		for (int i = 0; i < embeddingDimension; i++) {
			if (i >= syzygies.rank()) {
				result.add(base.zero());
			} else {
				result.add(diagonal.entry(i + 1, i + 1));
			}
		}
		return result;
	}

	public List<Ideal<T>> diagonalIdeals() {
		if (diagonalIdeals == null) {
			diagonalIdeals = new ArrayList<>();
			for (T diagonal : diagonalRanks()) {
				diagonalIdeals.add(base.getIdeal(Collections.singletonList(diagonal)));
			}
		}
		return diagonalIdeals;
	}

	private Vector<T> asDiagonalVector(Mod<S> s) {
		return asDiagonalVector(lift(s));
	}

	public Vector<T> asDiagonalVector(S s) {
		MatrixModule<T>.SmithNormalFormResult smith = matrixModule.smithNormalForm(syzygiesMatrix);
		return matrixModule.codomainAlgebra().multiply(smith.getRowOperationsInverse(), freeModule.asVector(s));
	}

	public Mod<S> fromDiagonalVector(Vector<T> t) {
		MatrixModule<T>.SmithNormalFormResult smith = matrixModule.smithNormalForm(syzygiesMatrix);
		List<T> reduced = new ArrayList<>();
		List<Ideal<T>> diagonalIdeals = diagonalIdeals();
		for (int i = 0; i < embeddingDimension; i++) {
			reduced.add(diagonalIdeals.get(i).residue(t.get(i + 1)));
		}
		return new Mod<>(freeModule
				.fromVector(matrixModule.codomainAlgebra().multiply(smith.getRowOperations(), new Vector<>(reduced))));
	}

//	private Mod<S> fromUnique(S s) {
//		if (rationalLattice != null) {
//			FiniteRationalVectorSpace space = rationalLattice.getVectorSpace();
//			Vector<T> asVector = freeModule.asVector(s);
//			Vector<Fraction> asRationalVector = Vector.mapVector(asFractionMap, asVector);
//			Vector<Fraction> closestSyzygy = space.closestLatticePoint(asRationalVector, rationalLattice);
//			Vector<Fraction> reduced = space.subtract(asRationalVector, closestSyzygy);
//			s = freeModule.fromVector(Vector.mapVector(fromFractionMap, reduced));
//		}
//		return new Mod<>(s);
//	}

	public List<Mod<S>> diagonalBasis() {
		FreeModule<T> freeModule = new FreeModule<>(base, embeddingDimension);
		List<Mod<S>> result = new ArrayList<>();
		for (Vector<T> basisVector : freeModule.getBasis()) {
			result.add(fromDiagonalVector(basisVector));
		}
		return result;
	}

	public Mod<S> reduce(S element) {
		Vector<T> inDiagonalBasis = asDiagonalVector(element);
		return fromDiagonalVector(inDiagonalBasis);
	}

	public S lift(Mod<S> element) {
		return element.representative;
	}
}
