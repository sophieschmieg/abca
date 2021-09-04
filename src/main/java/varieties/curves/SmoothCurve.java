package varieties.curves;

import java.util.List;

import fields.helper.CoordinateRing;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import varieties.ProjectivePoint;
import varieties.ProjectiveVariety;
import varieties.RationalFunction;
import varieties.Variety;
import varieties.curves.DivisorGroup.Divisor;

public interface SmoothCurve<T extends Element<T>> extends Variety<T, ProjectivePoint<T>> {
	public int getEmbeddingDimension();
	public ProjectiveVariety<T> asProjectiveVariety();
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p);
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p);
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div);
	public boolean isPrincipal(Divisor<T> div);
	public CoordinateRing<T> getCoordinateRing();
}
