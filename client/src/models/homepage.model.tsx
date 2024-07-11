import { IVUSSummary } from "./vus-summary.model";

export interface IHomepageData {
  vusList?: IVUSSummary[] | null;
  lastClinvarUpdateDate: string;
  lastPubUpdateDate: string;
}
