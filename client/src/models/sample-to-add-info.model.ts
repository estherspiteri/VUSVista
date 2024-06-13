import { IPhenotype } from "./phenotype.model";

export interface ISampleToAddInfo {
  sampleId: string;
  hgvs?: string | null;
  genotype?: string | null;
  phenotypes?: IPhenotype[] | null;
}
