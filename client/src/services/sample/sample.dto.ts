import { IPhenotype } from "../../models/phenotype.model";
import { ISample, ISampleSummary } from "../../models/view-samples.model";

export interface ILoadAllSamplesResponse {
  isSuccess: boolean;
  sampleList?: ISampleSummary[] | null;
}

export interface IGetSampleRequest {
  sampleId: string;
}

export interface IGetSampleResponse {
  isSuccess: boolean;
  sample: ISample;
}

export interface IGetHPOTermsRequest {
  phenotype: string;
}

export interface IHPOTerm {
  ontologyId: string;
  name: string;
}

export interface IGetHPOTermsResponse {
  isSuccess: boolean;
  hpoTerms: IHPOTerm[];
}

export interface IAddPhenotypeRequest {
  phenotype: IPhenotype;
  sampleId: string;
}

export interface IAddPhenotypeResponse {
  isSuccess: boolean;
}

export interface IRemovePhenotypeRequest {
  phenotype: IPhenotype;
  sampleId: string;
}

export interface IRemovePhenotypeResponse {
  isSuccess: boolean;
}
