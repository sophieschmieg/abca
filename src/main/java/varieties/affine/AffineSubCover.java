package varieties.affine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.interfaces.Element;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import varieties.affine.AffineMorphism.AffineFiberedProduct;

public class AffineSubCover<T extends Element<T>> {
	private AffineCover<T> coarse;
	private AffineCover<T> fine;
	private List<AffineMorphism<T>> refinement;
	private List<Integer> coarseCoverIndex;
	private List<List<Integer>> fineCoverIndex;

	public AffineSubCover(AffineCover<T> coarse, List<List<AffineMorphism<T>>> refinement) {
		this.coarse = coarse;
		this.refinement = new ArrayList<>();
		this.coarseCoverIndex = new ArrayList<>();
		this.fineCoverIndex = new ArrayList<>();
		List<AffineScheme<T>> fineCover = new ArrayList<>();
		List<List<AffineMorphism<T>>> fineAtlas = new ArrayList<>();
		for (int i = 0; i < coarse.getCover().size(); i++) {
			List<Integer> fineCoverIndeces = new ArrayList<>();
			CoordinateRing<T> coarseCoordinateRing = coarse.getCover().get(i).getCoordinateRing();
			CoordinateIdeal<T> coarseMissedIdeal = coarseCoordinateRing.getZeroIdeal();
			this.fineCoverIndex.add(fineCoverIndeces);
			for (int j = 0; j < refinement.get(i).size(); j++) {
				AffineMorphism<T> morphism = refinement.get(i).get(j);
				if (!morphism.getRange().equals(coarse.getCover().get(i))) {
					throw new ArithmeticException("Not a subcover!");
				}
				if (!morphism.isOpenImmersion()) {
					throw new ArithmeticException("Not an embedding!");
				}
				coarseMissedIdeal = (CoordinateIdeal<T>) coarseCoordinateRing.add(
						coarseCoordinateRing.getIdeal(Collections.singletonList(morphism.getInversionObstruction())),
						coarseMissedIdeal);
				int fineCoverIndex = fineCover.size();
				fineCoverIndeces.add(fineCoverIndex);
				fineCover.add(morphism.getDomain());
				this.refinement.add(morphism);
				this.coarseCoverIndex.add(i);
			}
			if (!coarseMissedIdeal.contains(coarseCoordinateRing.one())) {
				throw new ArithmeticException("Not a subcover, !");
			}
		}
		for (int i = 0; i < fineCover.size(); i++) {
			List<AffineMorphism<T>> intersections = new ArrayList<>();
			for (int j = 0; j < fineCover.size(); j++) {
				intersections.add(null);
			}
			intersections.set(i, fineCover.get(i).identityMorphism());
			fineAtlas.add(intersections);
		}
		for (int fineIndex = 1; fineIndex < fineCover.size(); fineIndex++) {
			int coarseIndex = coarseCoverIndex.get(fineIndex);
			AffineMorphism<T> morphism = this.refinement.get(fineIndex);
			for (int otherFineIndex = 0; otherFineIndex < fineIndex; otherFineIndex++) {
				int otherCoarseIndex = coarseCoverIndex.get(otherFineIndex);
				AffineMorphism<T> otherMorphism = this.refinement.get(otherFineIndex);
				AffineMorphism<T> coarseIntersection = coarse.getIntersectionEmbedding(coarseIndex, otherCoarseIndex);
				AffineFiberedProduct<T> refinementCoarseIntersection = AffineMorphism.fiberedProduct(coarseIntersection,
						morphism);
				AffineMorphism<T> otherCoarseIntersection = coarse.getIntersectionEmbedding(otherCoarseIndex,
						coarseIndex);
				AffineFiberedProduct<T> otherRefinementOtherCoarseIntersection = AffineMorphism
						.fiberedProduct(otherCoarseIntersection, otherMorphism);
				AffineFiberedProduct<T> refinedIntersection = AffineMorphism.fiberedProduct(
						refinementCoarseIntersection.getSimplifiedProjection1(),
						otherRefinementOtherCoarseIntersection.getSimplifiedProjection1());
				fineAtlas.get(fineIndex).set(otherFineIndex,
						AffineMorphism.composition(refinedIntersection.getSimplifiedProjection1(),
								refinementCoarseIntersection.getSimplifiedProjection2()));
				fineAtlas.get(otherFineIndex).set(fineIndex,
						AffineMorphism.composition(refinedIntersection.getSimplifiedProjection2(),
								otherRefinementOtherCoarseIntersection.getSimplifiedProjection2()));
			}
		}
		this.fine = new AffineCover<>(fineCover, fineAtlas);
	}

	private static <T extends Element<T>> List<List<AffineMorphism<T>>> trivialRefinement(AffineCover<T> cover) {
		List<List<AffineMorphism<T>>> result = new ArrayList<>();
		for (int i = 0; i < cover.getCover().size(); i++) {
			result.add(Collections.singletonList(cover.getCover().get(i).identityMorphism()));
		}
		return result;
	}

	public AffineSubCover(AffineCover<T> cover) {
		this(cover, trivialRefinement(cover));
	}

	public static <T extends Element<T>> AffineSubCover<T> composition(AffineSubCover<T> fine,
			AffineSubCover<T> coarse) {
		if (!fine.getCoarseCover().equals(coarse.getFineCover())) {
			throw new ArithmeticException("Not a refinement!");
		}
		List<List<AffineMorphism<T>>> refinement = new ArrayList<>();
		for (int i = 0; i < coarse.getCoarseCover().size(); i++) {
			List<AffineMorphism<T>> refine = new ArrayList<>();
			refinement.add(refine);
			for (int j : coarse.getFineCoverIndex(i)) {
				for (int k : fine.getFineCoverIndex(j)) {
					refine.add(AffineMorphism.composition(fine.getRefinement(k), coarse.getRefinement(j)));
				}
			}
		}
		return new AffineSubCover<>(coarse.getCoarseCover(), refinement);
	}

	public static class CommonRefinement<T extends Element<T>> {
		private AffineCover<T> cover;
		private AffineSubCover<T> firstSubCover;
		private AffineSubCover<T> secondSubCover;
		private AffineSubCover<T> refinement;
		private AffineSubCover<T> firstRefinement;
		private AffineSubCover<T> secondRefinement;

		public CommonRefinement(AffineSubCover<T> firstSubCover, AffineSubCover<T> secondSubCover) {
			if (!firstSubCover.getCoarseCover().equals(secondSubCover.getCoarseCover())) {
				throw new ArithmeticException("Different coarse cover");
			}
			this.cover = firstSubCover.getCoarseCover();
			this.firstSubCover = firstSubCover;
			this.secondSubCover = secondSubCover;
			List<List<AffineMorphism<T>>> commonRefinementList = new ArrayList<>();
			for (int i = 0; i < cover.size(); i++) {
				commonRefinementList.add(new ArrayList<>());
			}
			List<List<AffineMorphism<T>>> firstRefinementList = new ArrayList<>();
			for (int i = 0; i < firstSubCover.getFineCover().size(); i++) {
				firstRefinementList.add(new ArrayList<>());
			}
			List<List<AffineMorphism<T>>> secondRefinementList = new ArrayList<>();
			for (int i = 0; i < secondSubCover.getFineCover().size(); i++) {
				secondRefinementList.add(new ArrayList<>());
			}
			for (int i = 0; i < cover.size(); i++) {
				for (int j : firstSubCover.getFineCoverIndex(i)) {
					for (int k : secondSubCover.getFineCoverIndex(i)) {
						AffineFiberedProduct<T> fiberedProduct = AffineMorphism
								.fiberedProduct(firstSubCover.getRefinement(j), secondSubCover.getRefinement(k));
						commonRefinementList.get(i).add(AffineMorphism.composition(fiberedProduct.getProjection1(),
								fiberedProduct.getStructural1()));
						firstRefinementList.get(j).add(fiberedProduct.getProjection1());
						secondRefinementList.get(k).add(fiberedProduct.getProjection2());
					}
				}
			}

			this.refinement = new AffineSubCover<>(cover, commonRefinementList);
			this.firstRefinement = new AffineSubCover<>(firstSubCover.getFineCover(), firstRefinementList);
			this.secondRefinement = new AffineSubCover<>(secondSubCover.getFineCover(), secondRefinementList);
		}

		public AffineCover<T> getCover() {
			return cover;
		}

		public AffineSubCover<T> getFirstSubCover() {
			return firstSubCover;
		}

		public AffineSubCover<T> getSecondSubCover() {
			return secondSubCover;
		}

		public AffineSubCover<T> getRefinement() {
			return refinement;
		}

		public AffineSubCover<T> getFirstRefinement() {
			return firstRefinement;
		}

		public AffineSubCover<T> getSecondRefinement() {
			return secondRefinement;
		}
	}

	public AffineCover<T> getCoarseCover() {
		return coarse;
	}

	public AffineCover<T> getFineCover() {
		return fine;
	}

	public AffineMorphism<T> getRefinement(int fineCoverIndex) {
		return refinement.get(fineCoverIndex);
	}

	public int getCoarseCoverIndex(int fineCoverIndex) {
		return coarseCoverIndex.get(fineCoverIndex);
	}

	public List<Integer> getFineCoverIndex(int coarseCoverIndex) {
		return fineCoverIndex.get(coarseCoverIndex);
	}
}
