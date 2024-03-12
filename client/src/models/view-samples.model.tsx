import { IHPOTerm } from "../services/vus/vus.dto";

export enum Genotype {
  Homozygous = "HOMOZYGOUS",
  Heterozygous = "HETEROZYGOUS",
}

export interface ISampleVariant {
  variantId: number;
  genotype: Genotype;
  acmgRuleNames?: string[];
}

export interface ISample {
  sampleId: string;
  genomeVersion: string;
  fileUploadName: string;
  dateOfFileUpload: Date;
  variants: ISampleVariant[];
  phenotype?: IHPOTerm[] | null;
}

export interface ISampleSummary {
  sampleId: string;
  dateOfFileUpload: Date;
  numOfVariants: number;
}
