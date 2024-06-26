import customFetch from "../api/api.service";
import {
  IAddPhenotypeRequest,
  IAddPhenotypeResponse,
  IAddVariantsRequest,
  IAddVariantsResponse,
  IDeleteSampleRequest,
  IDeleteSampleResponse,
  IGetHPOTermsRequest,
  IGetHPOTermsResponse,
  IGetSampleRequest,
  IGetSampleResponse,
  ILoadAllSamplesResponse,
  IRemovePhenotypeRequest,
  IRemovePhenotypeResponse,
  IRemoveVariantsRequest,
  IRemoveVariantsResponse,
  IUpdateHgvsRequest,
  IUpdateHgvsResponse,
} from "./sample.dto";

export class SampleService {
  async loadAllSamples(): Promise<ILoadAllSamplesResponse> {
    const result: ILoadAllSamplesResponse = await customFetch(`/sample/view`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getSample(input: IGetSampleRequest): Promise<IGetSampleResponse> {
    const result: IGetSampleResponse = await customFetch(
      `/sample/view/${input.sampleId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getHPOTerms(input: IGetHPOTermsRequest): Promise<IGetHPOTermsResponse> {
    const result: IGetHPOTermsResponse = await customFetch(
      `/sample/phenotype/${input.phenotype}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addPhenotype(
    input: IAddPhenotypeRequest
  ): Promise<IAddPhenotypeResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("phenotype", JSON.stringify(input.phenotype));

    const result: IGetHPOTermsResponse = await customFetch(
      `/sample/add-phenotype`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async removePhenotype(
    input: IRemovePhenotypeRequest
  ): Promise<IRemovePhenotypeResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("phenotype", JSON.stringify(input.phenotype));

    const result: IGetHPOTermsResponse = await customFetch(
      `/sample/remove-phenotype`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async deleteSample(
    input: IDeleteSampleRequest
  ): Promise<IDeleteSampleResponse> {
    const result: IGetHPOTermsResponse = await customFetch(
      `/sample/delete/${input.sampleId}`,
      {
        method: "DELETE",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async updateHgvs(input: IUpdateHgvsRequest): Promise<IUpdateHgvsResponse> {
    const result: IUpdateHgvsResponse = await customFetch(
      `/sample/edit-hgvs/${input.sampleId}/${input.variantId}/${input.updatedHgvs}`,
      {
        method: "POST",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addVariants(input: IAddVariantsRequest): Promise<IAddVariantsResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("variantsToAdd", JSON.stringify(input.variantsToAdd));

    const result: IAddVariantsResponse = await customFetch(
      `/sample/add-variants/${input.sampleId}`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async removeVariants(
    input: IRemoveVariantsRequest
  ): Promise<IRemoveVariantsResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("deleteSample", JSON.stringify(input.isDeleteSample));
    data.append("variantIdsToRemove", JSON.stringify(input.variantIdsToRemove));

    const result: IRemoveVariantsResponse = await customFetch(
      `/sample/remove-variants/${input.sampleId}`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const samplesService = new SampleService();
