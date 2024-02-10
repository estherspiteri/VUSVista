import { IHPOTerm } from "../services/vus/vus.dto";

export enum Genotype {
  Homozygous = "HOMOZYGOUS",
  Heterozygous = "HETEROZYGOUS",
}

export interface ISampleVariant {
  variantId: number;
  genotype: Genotype;
}

export interface ISample {
  sampleId: string;
  genomeVersion: string;
  fileUploadName: string;
  dateOfFileUpload: Date;
  variants: ISampleVariant[];
  phenotype?: IHPOTerm[] | null;
}
