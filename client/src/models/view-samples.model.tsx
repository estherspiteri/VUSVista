import { IHPOTerm } from "../services/sample/sample.dto";
import { IVUSSummary } from "./publication-view.model";

export enum Genotype {
  Homozygous = "HOMOZYGOUS",
  Heterozygous = "HETEROZYGOUS",
}

export interface ISampleVariant {
  variantId: number;
  variant: IVUSSummary;
  genotype: Genotype;
  acmgRuleIds?: number[];
}

export interface IFile {
  filename: string;
  dateOfFileUpload: Date;
}

export interface ISample {
  sampleId: string;
  genomeVersion: string;
  files: IFile[];
  variants: ISampleVariant[];
  phenotype?: IHPOTerm[] | null;
}

export interface ISampleSummary {
  sampleId: string;
  numOfVariants: number;
}
