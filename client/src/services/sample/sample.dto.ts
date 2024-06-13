import { IPhenotype } from "../../models/phenotype.model";
import { IVariantToAddInfo } from "../../models/variant-to-add-info.model";
import {
  INotSampleVariant,
  ISample,
  ISampleSummary,
  ISampleVariant,
} from "../../models/view-samples.model";

export interface ILoadAllSamplesResponse {
  isSuccess: boolean;
  sampleList?: ISampleSummary[] | null;
}

export interface IGetSampleRequest {
  sampleId: string;
}

export interface IGetSampleResponse {
  isSuccess: boolean;
  sample?: ISample | null;
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

export interface IDeleteSampleRequest {
  sampleId: string;
}

export interface IDeleteSampleResponse {
  isSuccess: boolean;
}

export interface IUpdateHgvsRequest {
  sampleId: string;
  variantId: string;
  updatedHgvs: string;
}

export interface IUpdateHgvsResponse {
  isSuccess: boolean;
}

export interface IAddVariantsRequest {
  sampleId: string;
  variantsToAdd: IVariantToAddInfo[];
}

export interface IAddVariantsResponse {
  isSuccess: boolean;
  updatedVariants: ISampleVariant[];
  updatedNotSampleVariants: INotSampleVariant[];
}

export interface IRemoveVariantsRequest {
  sampleId: string;
  variantIdsToRemove: number[];
  isDeleteSample: boolean;
}

export interface IRemoveVariantsResponse {
  isSuccess: boolean;
  isSampleDeleted: boolean;
  updatedVariants?: ISampleVariant[] | null;
  updatedNotSampleVariants?: INotSampleVariant[] | null;
}
