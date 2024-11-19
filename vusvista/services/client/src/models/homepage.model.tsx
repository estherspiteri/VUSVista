import { IVUSSummary } from "./vus-summary.model";

export interface IHomepageData {
  vusList?: IVUSSummary[] | null;
  lastClinvarUpdateDate?: string | null;
  lastPubUpdateDate?: string | null;
}
