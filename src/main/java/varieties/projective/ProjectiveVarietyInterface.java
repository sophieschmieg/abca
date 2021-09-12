package varieties.projective;

import fields.interfaces.Element;
import varieties.Scheme;

public interface ProjectiveVarietyInterface<T extends Element<T>> extends Scheme<T, ProjectivePoint<T>> {
	ProjectiveScheme<T> asProjectiveVariety();
}
